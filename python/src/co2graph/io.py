"""Load Julia-exported graph data and snapshots."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class TrapNode:
    """Static properties of a single trap."""

    trap_id: int
    volume_capacity: float
    subvolume: float
    net_capacity: float
    spillpoint_elevation: float
    min_topography: float
    parent: int  # 0 = root (no parent)
    children: list[int]
    leakage_height: float  # Inf encoded as 1e308
    leakage_location: tuple[int, int]  # (row, col) 0-indexed, unpadded
    footprint_size: int


@dataclass
class LayerData:
    """All exported data for a single layer."""

    index: int
    name: str
    topography: np.ndarray  # (nx, ny)
    regions: np.ndarray  # (nx, ny) int, trap assignment per cell
    nodes: list[TrapNode]
    spill_targets: np.ndarray  # (n_traps,) parent trap per trap (0 = root)
    # Ragged footprints: cell (row, col) for each trap
    footprint_offsets: np.ndarray  # (n_traps + 1,)
    footprint_rows: np.ndarray  # flat
    footprint_cols: np.ndarray  # flat
    # Z-vol tables for volume-to-height conversion
    zvol_offsets: np.ndarray  # (n_traps + 1,)
    zvol_z: np.ndarray  # flat z values
    zvol_v: np.ndarray  # flat volume values

    @property
    def n_traps(self) -> int:
        return len(self.nodes)

    def get_footprint(self, trap_id: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (rows, cols) arrays for a trap's footprint (1-indexed trap_id)."""
        i = trap_id - 1
        s, e = self.footprint_offsets[i], self.footprint_offsets[i + 1]
        return self.footprint_rows[s:e], self.footprint_cols[s:e]

    def get_zvol_table(self, trap_id: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (z_values, vol_values) for a trap's z-vol table (1-indexed trap_id)."""
        i = trap_id - 1
        s, e = self.zvol_offsets[i], self.zvol_offsets[i + 1]
        return self.zvol_z[s:e], self.zvol_v[s:e]


@dataclass
class InterLayerOverlap:
    """Overlap between traps in adjacent layers."""

    src_layer: int
    dst_layer: int
    src_traps: np.ndarray  # trap IDs in lower layer
    dst_traps: np.ndarray  # trap IDs in upper layer
    overlap_counts: np.ndarray  # number of shared grid cells


@dataclass
class GraphData:
    """Complete exported graph data."""

    n_layers: int
    nx: int
    ny: int
    dx: float
    dy: float
    depth_min: float
    depth_max: float
    layer_names: list[str]
    boundary_condition: str
    pad_width: int
    reservoir_properties: list[dict]
    layers: list[LayerData]
    overlaps: list[InterLayerOverlap]


@dataclass
class SnapshotData:
    """Simulation snapshots for all layers."""

    timepoints: np.ndarray  # (T,)
    volumes: list[np.ndarray]  # per layer: (T, n_traps)
    filled: list[np.ndarray]  # per layer: (T, n_traps) bool
    leaking: list[np.ndarray]  # per layer: (T, n_traps) bool
    draining: list[np.ndarray]  # per layer: (T, n_traps) bool
    totals: list[np.ndarray]  # per layer: (T, 4) [injected, stored, drained, passthrough]


def load_graph_data(path: str | Path) -> GraphData:
    """Load graph data exported by Julia's ``export_graph_data``."""
    path = Path(path)

    with open(path / "metadata.json") as f:
        meta = json.load(f)

    n_layers = meta["n_layers"]
    layers: list[LayerData] = []

    for li in range(1, n_layers + 1):
        topo = np.load(path / f"layer_{li}_topography.npy")
        regions = np.load(path / f"layer_{li}_regions.npy")
        spill_targets = np.load(path / f"layer_{li}_spill_targets.npy")

        with open(path / f"layer_{li}_nodes.json") as f:
            raw_nodes = json.load(f)

        nodes = [
            TrapNode(
                trap_id=n["trap_id"],
                volume_capacity=n["volume_capacity"],
                subvolume=n["subvolume"],
                net_capacity=n["net_capacity"],
                spillpoint_elevation=n["spillpoint_elevation"],
                min_topography=n["min_topography"],
                parent=n["parent"],
                children=n["children"],
                leakage_height=n["leakage_height"] if n["leakage_height"] is not None else float("inf"),
                leakage_location=tuple(n["leakage_location"]),
                footprint_size=n["footprint_size"],
            )
            for n in raw_nodes
        ]

        fp_offsets = np.load(path / f"layer_{li}_footprint_offsets.npy")
        fp_rows = np.load(path / f"layer_{li}_footprint_rows.npy")
        fp_cols = np.load(path / f"layer_{li}_footprint_cols.npy")

        zv_offsets = np.load(path / f"layer_{li}_zvol_offsets.npy")
        zv_z = np.load(path / f"layer_{li}_zvol_z.npy")
        zv_v = np.load(path / f"layer_{li}_zvol_v.npy")

        layers.append(
            LayerData(
                index=li,
                name=meta["layer_names"][li - 1],
                topography=topo,
                regions=regions,
                nodes=nodes,
                spill_targets=spill_targets,
                footprint_offsets=fp_offsets,
                footprint_rows=fp_rows,
                footprint_cols=fp_cols,
                zvol_offsets=zv_offsets,
                zvol_z=zv_z,
                zvol_v=zv_v,
            )
        )

    # Load inter-layer overlaps
    overlaps: list[InterLayerOverlap] = []
    for li in range(1, n_layers):
        src_file = path / f"overlap_{li}_{li + 1}_src.npy"
        if src_file.exists():
            overlaps.append(
                InterLayerOverlap(
                    src_layer=li,
                    dst_layer=li + 1,
                    src_traps=np.load(src_file),
                    dst_traps=np.load(path / f"overlap_{li}_{li + 1}_dst.npy"),
                    overlap_counts=np.load(path / f"overlap_{li}_{li + 1}_count.npy"),
                )
            )

    return GraphData(
        n_layers=n_layers,
        nx=meta["nx"],
        ny=meta["ny"],
        dx=meta["dx"],
        dy=meta["dy"],
        depth_min=meta["depth_min"],
        depth_max=meta["depth_max"],
        layer_names=meta["layer_names"],
        boundary_condition=meta["boundary_condition"],
        pad_width=meta["pad_width"],
        reservoir_properties=meta["reservoir_properties"],
        layers=layers,
        overlaps=overlaps,
    )


def load_snapshots(path: str | Path, subdir: str = "snapshots") -> SnapshotData:
    """Load simulation snapshots exported by Julia's ``export_snapshots``."""
    path = Path(path) / subdir

    timepoints = np.load(path / "timepoints.npy")

    volumes, filled, leaking, draining, totals = [], [], [], [], []
    li = 1
    while (path / f"layer_{li}_volumes.npy").exists():
        volumes.append(np.load(path / f"layer_{li}_volumes.npy"))
        filled.append(np.load(path / f"layer_{li}_filled.npy"))
        leaking.append(np.load(path / f"layer_{li}_leaking.npy"))
        draining.append(np.load(path / f"layer_{li}_draining.npy"))
        totals.append(np.load(path / f"layer_{li}_totals.npy"))
        li += 1

    return SnapshotData(
        timepoints=timepoints,
        volumes=volumes,
        filled=filled,
        leaking=leaking,
        draining=draining,
        totals=totals,
    )
