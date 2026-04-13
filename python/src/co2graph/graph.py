"""Build GNN-compatible graphs from exported trap data."""

from __future__ import annotations

import numpy as np

from co2graph.io import GraphData, LayerData


def build_trap_graph(
    data: GraphData,
    layer_idx: int,
    *,
    node_feature_keys: list[str] | None = None,
) -> dict[str, np.ndarray]:
    """Build a directed graph on the trap hierarchy of a single layer.

    Edges follow spill direction: child -> parent (CO2 flows uphill in the
    hierarchy when a trap fills and spills). Each edge also has a reverse
    counterpart for bidirectional message passing.

    Parameters
    ----------
    data : GraphData
        Loaded export data.
    layer_idx : int
        1-indexed layer number.
    node_feature_keys : list[str] | None
        Which TrapNode fields to include in the node feature matrix.
        Defaults to a sensible set of geometric/capacity features.

    Returns
    -------
    dict with keys:
        ``senders``        — (E,) int32 source node indices (0-indexed)
        ``receivers``      — (E,) int32 target node indices (0-indexed)
        ``node_features``  — (N, F) float32 node feature matrix
        ``edge_features``  — (E, D) float32 edge feature matrix
        ``n_nodes``        — int
        ``n_edges``        — int
    """
    layer = data.layers[layer_idx - 1]
    n = layer.n_traps

    if node_feature_keys is None:
        node_feature_keys = [
            "volume_capacity",
            "net_capacity",
            "spillpoint_elevation",
            "min_topography",
            "leakage_height",
            "footprint_size",
        ]

    # ── Node features ──
    node_features = np.zeros((n, len(node_feature_keys)), dtype=np.float32)
    for i, node in enumerate(layer.nodes):
        for j, key in enumerate(node_feature_keys):
            val = float(getattr(node, key))
            # Clamp Inf to a large finite value for the GNN
            node_features[i, j] = min(val, 1e6) if np.isfinite(val) else 0.0

    # ── Edges: parent-child (bidirectional) ──
    senders_list: list[int] = []
    receivers_list: list[int] = []
    edge_feats_list: list[list[float]] = []

    for node in layer.nodes:
        t = node.trap_id - 1  # 0-indexed
        if node.parent > 0:
            p = node.parent - 1  # 0-indexed

            # Child → parent (spill direction)
            delta_elev = (
                node_features[p, 2] - node_features[t, 2]
            )  # spillpoint elevation diff
            senders_list.append(t)
            receivers_list.append(p)
            edge_feats_list.append([delta_elev, 1.0])  # type=1 (child→parent)

            # Parent → child (reverse)
            senders_list.append(p)
            receivers_list.append(t)
            edge_feats_list.append([-delta_elev, -1.0])  # type=-1 (parent→child)

    senders = np.array(senders_list, dtype=np.int32)
    receivers = np.array(receivers_list, dtype=np.int32)
    edge_features = (
        np.array(edge_feats_list, dtype=np.float32)
        if edge_feats_list
        else np.zeros((0, 2), dtype=np.float32)
    )

    return {
        "senders": senders,
        "receivers": receivers,
        "node_features": node_features,
        "edge_features": edge_features,
        "n_nodes": n,
        "n_edges": len(senders_list),
    }


def build_multi_layer_graph(
    data: GraphData,
    *,
    node_feature_keys: list[str] | None = None,
    include_leakage_edges: bool = True,
) -> dict[str, np.ndarray]:
    """Build a single graph spanning all layers with inter-layer leakage edges.

    Nodes are traps across all layers. Intra-layer edges come from spill
    hierarchy; inter-layer edges connect overlapping traps (leakage paths).

    Node indices are globally offset: layer 1 traps are 0..N1-1,
    layer 2 traps are N1..N1+N2-1, etc.

    Parameters
    ----------
    data : GraphData
        Loaded export data.
    node_feature_keys : list[str] | None
        Feature keys (see ``build_trap_graph``).
    include_leakage_edges : bool
        Whether to add edges between overlapping traps in adjacent layers.

    Returns
    -------
    dict with keys: ``senders``, ``receivers``, ``node_features``,
    ``edge_features``, ``n_nodes``, ``n_edges``, ``layer_offsets``,
    ``layer_ids``.
    """
    if node_feature_keys is None:
        node_feature_keys = [
            "volume_capacity",
            "net_capacity",
            "spillpoint_elevation",
            "min_topography",
            "leakage_height",
            "footprint_size",
        ]

    # Build per-layer graphs and accumulate
    all_node_features: list[np.ndarray] = []
    all_senders: list[np.ndarray] = []
    all_receivers: list[np.ndarray] = []
    all_edge_features: list[np.ndarray] = []
    layer_offsets = [0]
    layer_ids: list[np.ndarray] = []

    for li in range(1, data.n_layers + 1):
        g = build_trap_graph(data, li, node_feature_keys=node_feature_keys)
        offset = layer_offsets[-1]

        all_node_features.append(g["node_features"])
        all_senders.append(g["senders"] + offset)
        all_receivers.append(g["receivers"] + offset)
        all_edge_features.append(g["edge_features"])
        layer_ids.append(np.full(g["n_nodes"], li, dtype=np.int32))
        layer_offsets.append(offset + g["n_nodes"])

    # Inter-layer leakage edges
    if include_leakage_edges:
        for ov in data.overlaps:
            offset_src = layer_offsets[ov.src_layer - 1]
            offset_dst = layer_offsets[ov.dst_layer - 1]

            # Leakage goes upward: lower layer trap → upper layer trap (bidirectional)
            src_global = ov.src_traps - 1 + offset_src  # 0-indexed
            dst_global = ov.dst_traps - 1 + offset_dst

            n_edges = len(src_global)
            overlap_frac = ov.overlap_counts.astype(np.float32)

            # Lower → upper (leakage direction)
            all_senders.append(src_global.astype(np.int32))
            all_receivers.append(dst_global.astype(np.int32))
            all_edge_features.append(
                np.column_stack([overlap_frac, np.full(n_edges, 2.0, dtype=np.float32)])
            )

            # Upper → lower (reverse)
            all_senders.append(dst_global.astype(np.int32))
            all_receivers.append(src_global.astype(np.int32))
            all_edge_features.append(
                np.column_stack(
                    [overlap_frac, np.full(n_edges, -2.0, dtype=np.float32)]
                )
            )

    return {
        "senders": np.concatenate(all_senders)
        if all_senders
        else np.array([], dtype=np.int32),
        "receivers": np.concatenate(all_receivers)
        if all_receivers
        else np.array([], dtype=np.int32),
        "node_features": np.concatenate(all_node_features, axis=0),
        "edge_features": np.concatenate(all_edge_features, axis=0)
        if all_edge_features
        else np.zeros((0, 2), dtype=np.float32),
        "n_nodes": layer_offsets[-1],
        "n_edges": sum(len(s) for s in all_senders),
        "layer_offsets": np.array(layer_offsets, dtype=np.int32),
        "layer_ids": np.concatenate(layer_ids),
    }


def trap_volumes_to_grid(
    layer: LayerData,
    volumes: np.ndarray,
    *,
    leaking: np.ndarray | None = None,
) -> np.ndarray:
    """Convert per-trap volumes to a 2D CO2 column-height grid.

    This replicates the Julia ``height_map`` logic so that GNN predictions
    (which output per-trap volumes) can be visualized on the spatial grid.

    Parameters
    ----------
    layer : LayerData
        Layer with topography, footprints, and z-vol tables.
    volumes : np.ndarray
        (n_traps,) volume per trap.
    leaking : np.ndarray | None
        (n_traps,) bool — if True, height is clamped to leakage_height.

    Returns
    -------
    np.ndarray
        (nx, ny) CO2 column height grid.
    """
    nx, ny = layer.topography.shape
    height_grid = np.zeros((nx, ny), dtype=np.float64)

    for i, node in enumerate(layer.nodes):
        vol = volumes[i]
        if vol <= 0.0:
            continue

        trap_id = node.trap_id

        # Determine CO2 column height
        if leaking is not None and leaking[i] and np.isfinite(node.leakage_height):
            h = node.leakage_height
        else:
            h = _volume_to_height(layer, trap_id, vol)

        if h <= 0.0:
            continue

        water_level = node.min_topography + h
        rows, cols = layer.get_footprint(trap_id)
        for r, c in zip(rows, cols):
            cell_h = max(0.0, water_level - layer.topography[r, c])
            height_grid[r, c] = max(height_grid[r, c], cell_h)

    return height_grid


def _volume_to_height(layer: LayerData, trap_id: int, volume: float) -> float:
    """Interpolate volume → CO2 column height using the z-vol table."""
    z_vals, v_vals = layer.get_zvol_table(trap_id)

    if len(z_vals) <= 1 or volume <= 0.0:
        return 0.0

    # Linear interpolation: volume → z (water level), then z - min_topo = height
    water_level = np.interp(volume, v_vals, z_vals)
    min_topo = layer.nodes[trap_id - 1].min_topography
    return max(0.0, water_level - min_topo)
