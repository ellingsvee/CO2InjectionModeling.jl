using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, parentof, subtrapsof,
    _compute_z_vol_tables

export export_graph_data, export_snapshots, build_graph

"""
    export_graph_data(layers, domain, reservoir_properties; path="graph_export")

Export the trap graph structure and topographic data for all layers to a directory
that can be loaded by the Python `co2graph` package.

Requires `NPZ` and `JSON` packages (e.g. `using NPZ, JSON`).

Creates:
- `metadata.json` — grid dimensions, layer info, reservoir properties
- `layer_{i}_topography.npy` — (nx, ny) unpadded topography
- `layer_{i}_regions.npy` — (nx, ny) unpadded trap-region map (0 = no trap)
- `layer_{i}_nodes.json` — per-trap features (capacity, spillpoint, parent, etc.)
- `layer_{i}_spill_targets.npy` — (n_traps,) spill target per trap (0 = out-of-domain)
- `layer_{i}_footprint_offsets.npy` + `layer_{i}_footprint_rows/cols.npy` — ragged footprints
- `layer_{i}_zvol_offsets.npy` + `layer_{i}_zvol_z/v.npy` — z-vol tables
- `overlap_{i}_{i+1}_*.npy` — inter-layer trap overlaps

# Arguments
- `layers`: `Vector{Layer}` from `analyze_base_surfaces`
- `domain`: `Domain3D` from `create_domain`
- `reservoir_properties`: `Vector{ReservoirProperties}` per layer (deepest first)
- `path`: output directory (created if needed)
"""
function export_graph_data end


"""
    export_snapshots(snapshots; path="graph_export", subdir="snapshots")

Export simulation snapshots (from `generate_multi_layer_snapshots`) to NPY files
for loading in Python as GNN training data.

Requires `NPZ` and `JSON` packages (e.g. `using NPZ, JSON`).

Creates `{path}/{subdir}/` with:
- `timepoints.npy` — (T,) vector of timestamps
- `layer_{i}_volumes.npy` — (T, n_traps) trap volumes per timepoint
- `layer_{i}_filled.npy` — (T, n_traps) bool
- `layer_{i}_leaking.npy` — (T, n_traps) bool
- `layer_{i}_draining.npy` — (T, n_traps) bool
- `layer_{i}_totals.npy` — (T, 4) [injected, stored, drained, passthrough]
"""
function export_snapshots end


"""
    build_graph(layers, reservoir_properties; include_leakage_edges=true) -> NamedTuple

Build a GNN-compatible graph spanning all layers.

Returns a `NamedTuple` with:
- `senders::Vector{Int32}` — 0-indexed source node indices
- `receivers::Vector{Int32}` — 0-indexed target node indices
- `node_features::Matrix{Float32}` — (n_nodes, n_features) feature matrix
- `edge_features::Matrix{Float32}` — (n_edges, n_edge_features) feature matrix
- `n_nodes::Int` — total number of nodes
- `n_edges::Int` — total number of edges
- `layer_offsets::Vector{Int32}` — start index of each layer's nodes (0-indexed)
- `layer_ids::Vector{Int32}` — layer ID per node (1-indexed)

Node features per trap (columns):
1. `volume_capacity` — total trap volume
2. `net_capacity` — capacity minus child subvolumes
3. `spillpoint_elevation` — water level at which trap spills
4. `min_topography` — lowest point in footprint
5. `leakage_height` — CO2 column threshold for leakage (clamped to 1e6)
6. `footprint_size` — number of grid cells in footprint

Edge features (columns):
1. `delta_elevation` — spillpoint elevation difference (sender→receiver)
2. `edge_type` — +1 child→parent, -1 parent→child, +2 lower→upper layer, -2 upper→lower

# Arguments
- `layers`: `Vector{Layer}` from `analyze_base_surfaces`
- `reservoir_properties`: `Vector{ReservoirProperties}` per layer (deepest first)
- `include_leakage_edges`: include inter-layer edges from trap overlap (default `true`)
"""
function build_graph(
    layers::Vector{Layer},
    reservoir_properties::Vector{ReservoirProperties};
    include_leakage_edges::Bool=true,
)
    n_layers = length(layers)

    # Collect per-layer node features and intra-layer edges
    all_node_features = Vector{Matrix{Float32}}()
    all_senders = Vector{Int32}()
    all_receivers = Vector{Int32}()
    all_edge_features = Vector{Vector{Float32}}()
    layer_offsets = Int32[0]
    layer_ids = Int32[]

    for (li, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        n_traps = numtraps(tstruct)
        rp = reservoir_properties[li]
        lh = rp.leakage_height
        offset = layer_offsets[end]

        # Node features: (n_traps, 6)
        nf = zeros(Float32, n_traps, 6)
        for t in 1:n_traps
            nf[t, 1] = Float32(tstruct.trapvolumes[t])
            nf[t, 2] = Float32(tstruct.trapvolumes[t] - tstruct.subvolumes[t])
            nf[t, 3] = Float32(tstruct.spillpoints[t].elevation)
            nf[t, 4] = Float32(get_min_topography_elevation(t, tstruct))
            raw_lh = isa(lh, Float64) ? lh : lh[t]
            nf[t, 5] = Float32(isfinite(raw_lh) ? raw_lh : Float32(1e6))
            nf[t, 6] = Float32(length(tstruct.footprints[t]))
        end
        push!(all_node_features, nf)

        # Intra-layer edges: bidirectional parent ↔ child
        for t in 1:n_traps
            p = parentof(tstruct, t)
            isnothing(p) && continue
            ti = Int32(t - 1 + offset)   # 0-indexed
            pi = Int32(p - 1 + offset)
            delta = nf[p, 3] - nf[t, 3]  # spillpoint elevation diff

            # child → parent
            push!(all_senders, ti)
            push!(all_receivers, pi)
            push!(all_edge_features, Float32[delta, 1.0f0])

            # parent → child
            push!(all_senders, pi)
            push!(all_receivers, ti)
            push!(all_edge_features, Float32[-delta, -1.0f0])
        end

        append!(layer_ids, fill(Int32(li), n_traps))
        push!(layer_offsets, Int32(offset + n_traps))
    end

    # Inter-layer leakage edges from trap overlap
    if include_leakage_edges
        nx, ny = get_grid_size(layers)
        for li in 1:n_layers-1
            pad_lo = layers[li].pad_width
            pad_hi = layers[li+1].pad_width
            regions_lo = layers[li].trap_structure.regions[pad_lo+1:end-pad_lo, pad_lo+1:end-pad_lo]
            regions_hi = layers[li+1].trap_structure.regions[pad_hi+1:end-pad_hi, pad_hi+1:end-pad_hi]

            offset_lo = layer_offsets[li]
            offset_hi = layer_offsets[li+1]

            # Find overlapping trap pairs
            overlap = Dict{Tuple{Int32,Int32}, Int32}()
            for idx in CartesianIndices(regions_lo)
                t_lo = regions_lo[idx]
                t_hi = regions_hi[idx]
                (t_lo > 0 && t_hi > 0) || continue
                key = (Int32(t_lo), Int32(t_hi))
                overlap[key] = get(overlap, key, Int32(0)) + Int32(1)
            end

            for ((t_lo, t_hi), count) in overlap
                si = Int32(t_lo - 1 + offset_lo)
                ri = Int32(t_hi - 1 + offset_hi)
                fcount = Float32(count)

                # lower → upper (leakage direction)
                push!(all_senders, si)
                push!(all_receivers, ri)
                push!(all_edge_features, Float32[fcount, 2.0f0])

                # upper → lower (reverse)
                push!(all_senders, ri)
                push!(all_receivers, si)
                push!(all_edge_features, Float32[fcount, -2.0f0])
            end
        end
    end

    # Assemble matrices
    node_features = vcat(all_node_features...)
    n_edges = length(all_senders)
    edge_features = n_edges > 0 ? reduce(hcat, all_edge_features)' : zeros(Float32, 0, 2)
    # Convert to a proper Matrix (reduce(hcat,...) gives an Adjoint)
    edge_features = Matrix{Float32}(edge_features)

    return (
        senders=all_senders,
        receivers=all_receivers,
        node_features=node_features,
        edge_features=edge_features,
        n_nodes=size(node_features, 1),
        n_edges=n_edges,
        layer_offsets=layer_offsets,
        layer_ids=layer_ids,
    )
end
