module CO2BatchFillExportExt

using CO2BatchFill
using NPZ
using JSON
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, parentof, subtrapsof,
    _compute_z_vol_tables

import CO2BatchFill: export_graph_data, export_snapshots


"""
    _write_npy(path, array)

Write an array to `.npy` format via NPZ.  NPZ.jl handles the column-major →
row-major conversion internally, so the shape in NumPy matches `size(array)`.
"""
function _write_npy(path::String, arr::AbstractArray)
    npzwrite(path, arr)
end


function export_graph_data(
    layers::Vector{CO2BatchFill.Layer},
    domain::CO2BatchFill.Domain3D,
    reservoir_properties::Vector{CO2BatchFill.ReservoirProperties};
    path::String="graph_export",
)
    mkpath(path)
    n_layers = length(layers)
    nx, ny = get_grid_size(layers)

    # ── metadata.json ──
    metadata = Dict(
        "n_layers" => n_layers,
        "nx" => nx,
        "ny" => ny,
        "dx" => domain.dx,
        "dy" => domain.dy,
        "depth_min" => domain.depth_min,
        "depth_max" => domain.depth_max,
        "layer_names" => [l.name for l in layers],
        "boundary_condition" => string(layers[1].boundary_condition),
        "pad_width" => layers[1].pad_width,
        "reservoir_properties" => [
            Dict(
                "sand_porosity" => rp.sand_porosity,
                "sand_residual_co2_saturation" => rp.sand_residual_co2_saturation,
                "sand_irreducible_water_saturation" => rp.sand_irreducible_water_saturation,
                "shale_pressure_threshold" => rp.shale_pressure_threshold,
                "leakage_height" => rp.leakage_height,
                "residual_leakage_time" => rp.residual_leakage_time,
                "brine_density" => rp.brine_density,
                "co2_density" => rp.co2_density,
            ) for rp in reservoir_properties
        ],
    )
    open(joinpath(path, "metadata.json"), "w") do io
        JSON.print(io, metadata, 2)
    end

    # ── Per-layer data ──
    for (li, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        pad = layer.pad_width
        n_traps = numtraps(tstruct)
        z_vol_tables = _compute_z_vol_tables(tstruct)

        # Topography (unpadded)
        topo = Float64.(tstruct.topography[pad+1:end-pad, pad+1:end-pad])
        _write_npy(joinpath(path, "layer_$(li)_topography.npy"), topo)

        # Regions (unpadded)
        regions = Int32.(tstruct.regions[pad+1:end-pad, pad+1:end-pad])
        _write_npy(joinpath(path, "layer_$(li)_regions.npy"), regions)

        # Spill targets: parent trap for each trap (0 = root)
        spill_targets = Int32[something(parentof(tstruct, t), 0) for t in 1:n_traps]
        _write_npy(joinpath(path, "layer_$(li)_spill_targets.npy"), spill_targets)

        # Node features as JSON
        rp = reservoir_properties[li]
        lh = rp.leakage_height
        cart_indices = CartesianIndices(tstruct.topography)
        nodes = []
        for t in 1:n_traps
            children = subtrapsof(tstruct, t)
            min_topo = get_min_topography_elevation(t, tstruct)
            sp_elev = tstruct.spillpoints[t].elevation

            leak_loc = find_leakage_location(t, tstruct)
            leak_i = leak_loc[1] - pad
            leak_j = leak_loc[2] - pad

            push!(nodes, Dict(
                "trap_id" => t,
                "volume_capacity" => tstruct.trapvolumes[t],
                "subvolume" => tstruct.subvolumes[t],
                "net_capacity" => tstruct.trapvolumes[t] - tstruct.subvolumes[t],
                "spillpoint_elevation" => sp_elev,
                "min_topography" => min_topo,
                "parent" => something(parentof(tstruct, t), 0),
                "children" => collect(children),
                "leakage_height" => isa(lh, Float64) ? lh : lh[t],
                "leakage_location" => [leak_i, leak_j],
                "footprint_size" => length(tstruct.footprints[t]),
            ))
        end
        open(joinpath(path, "layer_$(li)_nodes.json"), "w") do io
            JSON.print(io, nodes, 2)
        end

        # Footprints as ragged array (offsets + flat row/col indices)
        # Footprints are linear indices into the padded topography matrix.
        offsets = Int64[0]
        all_rows = Int32[]
        all_cols = Int32[]
        for t in 1:n_traps
            for lin_idx in tstruct.footprints[t]
                ci = cart_indices[lin_idx]
                push!(all_rows, Int32(ci[1] - 1 - pad))  # 0-indexed, unpadded
                push!(all_cols, Int32(ci[2] - 1 - pad))
            end
            push!(offsets, length(all_rows))
        end
        _write_npy(joinpath(path, "layer_$(li)_footprint_offsets.npy"), offsets)
        _write_npy(joinpath(path, "layer_$(li)_footprint_rows.npy"), all_rows)
        _write_npy(joinpath(path, "layer_$(li)_footprint_cols.npy"), all_cols)

        # Z-vol tables as ragged array
        zvol_offsets = Int64[0]
        all_z = Float64[]
        all_v = Float64[]
        for t in 1:n_traps
            zvals, vvals = z_vol_tables[t]
            append!(all_z, zvals)
            append!(all_v, vvals)
            push!(zvol_offsets, length(all_z))
        end
        _write_npy(joinpath(path, "layer_$(li)_zvol_offsets.npy"), zvol_offsets)
        _write_npy(joinpath(path, "layer_$(li)_zvol_z.npy"), all_z)
        _write_npy(joinpath(path, "layer_$(li)_zvol_v.npy"), all_v)
    end

    # Inter-layer overlap
    for li in 1:n_layers-1
        pad_lo = layers[li].pad_width
        pad_hi = layers[li+1].pad_width
        regions_lo = layers[li].trap_structure.regions[pad_lo+1:end-pad_lo, pad_lo+1:end-pad_lo]
        regions_hi = layers[li+1].trap_structure.regions[pad_hi+1:end-pad_hi, pad_hi+1:end-pad_hi]

        n_lo = numtraps(layers[li].trap_structure)

        overlap_src = Int32[]
        overlap_dst = Int32[]
        overlap_count = Int32[]

        for t_lo in 1:n_lo
            counts = Dict{Int32,Int32}()
            for idx in CartesianIndices(regions_lo)
                if regions_lo[idx] == t_lo && regions_hi[idx] > 0
                    t_hi = Int32(regions_hi[idx])
                    counts[t_hi] = get(counts, t_hi, Int32(0)) + Int32(1)
                end
            end
            for (t_hi, c) in counts
                push!(overlap_src, Int32(t_lo))
                push!(overlap_dst, Int32(t_hi))
                push!(overlap_count, c)
            end
        end
        _write_npy(joinpath(path, "overlap_$(li)_$(li+1)_src.npy"), overlap_src)
        _write_npy(joinpath(path, "overlap_$(li)_$(li+1)_dst.npy"), overlap_dst)
        _write_npy(joinpath(path, "overlap_$(li)_$(li+1)_count.npy"), overlap_count)
    end

    println("Exported graph data to: $path")
    println("  Layers: $n_layers, Grid: $(nx)×$(ny)")
    for (li, layer) in enumerate(layers)
        println("  Layer $li ($(layer.name)): $(numtraps(layer.trap_structure)) traps")
    end
    return path
end


function export_snapshots(
    snapshots::Vector{CO2BatchFill.MultiLayerSnapshot};
    path::String="graph_export",
    subdir::String="snapshots",
)
    outdir = joinpath(path, subdir)
    mkpath(outdir)

    T = length(snapshots)
    n_layers = length(snapshots[1].layers)

    timepoints = Float64[s.layers[1].timestamp for s in snapshots]
    _write_npy(joinpath(outdir, "timepoints.npy"), timepoints)

    for li in 1:n_layers
        n_traps = length(snapshots[1].layers[li].trap_volumes)

        volumes = zeros(Float64, T, n_traps)
        filled = zeros(Bool, T, n_traps)
        leaking = zeros(Bool, T, n_traps)
        draining = zeros(Bool, T, n_traps)
        totals = zeros(Float64, T, 4)

        for (ti, snap) in enumerate(snapshots)
            ls = snap.layers[li]
            volumes[ti, :] .= ls.trap_volumes
            filled[ti, :] .= ls.trap_filled
            leaking[ti, :] .= ls.trap_leaking
            draining[ti, :] .= ls.trap_draining
            totals[ti, :] .= [ls.total_injected, ls.total_stored, ls.total_drained, ls.total_passthrough]
        end

        _write_npy(joinpath(outdir, "layer_$(li)_volumes.npy"), volumes)
        _write_npy(joinpath(outdir, "layer_$(li)_filled.npy"), filled)
        _write_npy(joinpath(outdir, "layer_$(li)_leaking.npy"), leaking)
        _write_npy(joinpath(outdir, "layer_$(li)_draining.npy"), draining)
        _write_npy(joinpath(outdir, "layer_$(li)_totals.npy"), totals)
    end

    println("Exported $(T) snapshots for $(n_layers) layers to: $outdir")
    return outdir
end

end
