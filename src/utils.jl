import Graphs

export reconstruct_3d_lithology, create_trap_mask_3d, scale_unit_volume_to_physical, get_all_parents

"""
Reconstruct 3D lithology grid from topography surfaces.

Returns a 3D array where:
- 0 = caprock (impermeable)
- 1 = sand (permeable reservoir)
- 2 = shale (low permeability)

Parameters:
- topography: SleipnerTopography struct
- domain: Domain3D struct
"""
function reconstruct_3d_lithology(topography::SleipnerTopography, domain::Domain3D)
    nx, ny, nz = domain.nx, domain.ny, domain.nz
    dz = domain.dz
    depth_min = domain.depth_min
    depth_max = domain.depth_max

    println("\nReconstructing 3D lithology grid...")
    println("  Grid dimensions: $(nx) × $(ny) × $(nz)")
    println("  Depth range: $(depth_min) to $(depth_max) m")
    println("  Cell size: dx=$(domain.dx) m, dy=$(domain.dy) m, dz=$(dz) m")

    # Initialize with shale (default)
    lithology = fill(2, nx, ny, nz)

    # Create depth array for all k indices (vectorized)
    # k=1 is deepest, k=nz is shallowest
    depth_at_k = depth_max .- (0.5:nz) .* dz  # Vector of length nz

    # Process caprock layer (vectorized over k)
    println("  Processing caprock...")
    caprock_count = 0
    for k in 1:nz
        cell_depth = depth_at_k[k]
        # topography.top_caprock is already (nx, ny) - no transpose needed
        caprock_mask = cell_depth .< topography.top_caprock  # Broadcasting, result is (nx, ny)
        n_caprock_at_k = count(caprock_mask)
        caprock_count += n_caprock_at_k
        if n_caprock_at_k > 0
            lithology[caprock_mask, k] .= 0  # Caprock
        end
    end
    println("    Assigned $(caprock_count) caprock cells")

    # Process sand layers (vectorized)
    println("  Processing sand layers...")
    sand_count = 0
    for layer in topography.sand_layers
        # Layer surfaces are already (nx, ny) - no transpose needed
        layer_top = layer["top"]
        layer_base = layer["base"]

        layer_sand_count = 0
        for k in 1:nz
            cell_depth = depth_at_k[k]
            # Mask for cells in this sand layer at this depth
            in_layer = (layer_top .<= cell_depth) .& (cell_depth .<= layer_base)

            n_in_layer = count(in_layer)
            if n_in_layer > 0
                lithology[in_layer, k] .= 1  # Sand
                layer_sand_count += n_in_layer
            end
        end
        sand_count += layer_sand_count
    end
    println("    Assigned $(sand_count) sand cells")

    total_shale = prod(size(lithology)) - caprock_count - sand_count
    println("    Remaining $(total_shale) shale cells")

    return lithology
end

"""
Create a 3D mask for a specific trap within a layer.

This function maps a 2D trap footprint from the TrapStructure to the full 3D domain.
The mask will be true for all cells that belong to the trap in the vertical extent of the layer.

Parameters:
- layer: Layer struct containing trap_structure
- trap_id: The trap ID (index) to create mask for
- topography: SleipnerTopography struct
- layer_dict: Dictionary entry from topography.sand_layers for this layer
- dz: vertical grid spacing (meters)

Returns:
- mask_3d: 3D boolean array (nx, ny, nz) where true indicates cells belonging to the trap
"""
function create_trap_mask_3d(
    layer::Layer,
    trap_id::Int,
    domain::Domain3D,
)
    nx, ny = domain.nx, domain.ny
    dz = domain.dz
    depth_max = domain.depth_max
    nz = domain.nz

    topography_2d = layer.trap_structure.topography

    # Get the footprint of the trap (linear indices in 2D)
    trap_footprint = layer.trap_structure.footprints[trap_id]
    trap_footprint_2d = CartesianIndices(topography_2d)[trap_footprint]

    # Spillpoint elevation is the TOP (shallowest depth) where trap spills
    trap_bottom_elevation = layer.trap_structure.spillpoints[trap_id].elevation

    # Create depth array for all k indices (k=1 is deepest, k=nz is shallowest)
    depths = depth_max .- (0.5:nz) .* dz
    depths_3d = reshape(depths, 1, 1, nz)

    # Create height mask: cells between topography base and spillpoint elevation
    height_mask_3d = (depths_3d .<= trap_bottom_elevation) .& (depths_3d .>= reshape(topography_2d, nx, ny, 1))

    # Create 2D mask for the trap footprint
    footprint_mask_2d = falses(nx, ny)
    for idx in trap_footprint_2d
        footprint_mask_2d[idx.I[1], idx.I[2]] = true
    end

    # Broadcast footprint to 3D and combine with height mask
    footprint_mask_3d = reshape(footprint_mask_2d, nx, ny, 1)
    final_mask_3d = footprint_mask_3d .& height_mask_3d

    return final_mask_3d
end

"""
Create a 3D mask for multiple traps within a layer.

This function maps 2D trap footprints from the TrapStructure to the full 3D domain.
The mask will be true for all cells that belong to any of the specified traps.

Parameters:
- layer: Layer struct containing trap_structure
- trap_ids: Vector of trap IDs (indices) to create mask for
- domain: Domain3D struct

Returns:
- mask_3d: 3D boolean array (nx, ny, nz) where true indicates cells belonging to any of the traps
"""
function create_trap_mask_3d(
    layer::Layer,
    trap_ids::Vector{Int},
    domain::Domain3D,
)
    nx, ny = domain.nx, domain.ny
    nz = domain.nz

    # Initialize the combined mask
    combined_mask = falses(nx, ny, nz)

    # Add each trap to the mask
    for trap_id in trap_ids
        trap_mask = create_trap_mask_3d(layer, trap_id, domain)
        combined_mask .|= trap_mask
    end

    return combined_mask
end

# function simulation_layer_snapshots_from_spill_events(seq::Vector{SpillEvent}, timepoints::Vector{Float64}, tstruct::TrapStructure, reservoir_properties::ReservoirProperties, domain::Domain3D)
function simulation_layer_snapshots_from_spill_events(
    layer::Layer,
    seq::Vector{SpillEvent},
    domain::Domain3D,
    reservoir_properties::ReservoirProperties,
    injection_events::Vector{InjectionEvent};
    num_snapshots::Int,
    start_time::Float64,
    end_time::Float64
)
    tstruct = layer.trap_structure
    timepoints = collect(range(start_time, stop=end_time, length=num_snapshots))

    tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)
    water_content = [e[2] for e in tstates]

    total_contents = zeros(Float64, length(timepoints))
    for time_ix = 1:length(timepoints)
        for trap_ix = 1:numtraps(tstruct)
            content = water_content[time_ix][trap_ix]
            total_contents[time_ix] += content
        end
    end
    
    # TODO: Compute residual_trapped_co2_volume properly
    # Boolean vector indicating which traps have residual trapping
    # Length the same as the spill_event.filled, but all false for now
    residual_trapped_co2_volume = 0.0
    residual_trapped = Vector{Bool}(falses(numtraps(tstruct)))
    snapshots = Vector{SimulationLayerSnapshot}()
    seq_ix = 1  # Current position in the spill event sequence


    # Compute the height matrices for each timepoint
    height_matrices = compute_co2_height_matrix(
        seq,
        tstruct;
        timepoint = timepoints,
        use_layer_base = true,
        tstates = tstates
    )


    for time_ix = 1:length(timepoints)
        tp = timepoints[time_ix]

        # Find the correct sequence index for this timepoint
        # Advance seq_ix until we find the last event at or before this timepoint
        while seq_ix < length(seq) && seq[seq_ix + 1].timestamp <= tp
            seq_ix += 1
        end

        # Verify we found a valid position
        @assert seq[seq_ix].timestamp <= tp "Timepoint $(tp) is before first sequence event at $(seq[seq_ix].timestamp)"

        # Use the spill event at the correct sequence position
        spill_event = seq[seq_ix]

        total_co2_volume = swim_volume_to_physical_volume(total_contents[time_ix], reservoir_properties, domain)

        injected_volume = compute_total_injected_amount(injection_events, tp)


        push!(snapshots, SimulationLayerSnapshot(
            timepoints[time_ix],
            spill_event,
            height_matrices[time_ix],
            tstates[time_ix][1],
            injected_volume,
            total_co2_volume,
            total_co2_volume,  # Placeholder: all CO2 is mobile for now
            residual_trapped_co2_volume, # Placeholder
            residual_trapped # Placeholder
        ));

    end

    return snapshots
end

function compute_total_injected_amount(injection_events::Vector{InjectionEvent}, time::Float64)
    total = 0.0
    for event_idx in 1:length(injection_events)
        event = injection_events[event_idx]
        t_start = event.timestamp
        t_end = event_idx < length(injection_events) ? injection_events[event_idx + 1].timestamp : time

        # If the interval is entirely after the requested time, skip
        if time < t_start
            continue
        end

        # Only integrate up to 'time'
        interval_end = min(t_end, time)
        dt = interval_end - t_start
        if dt > 0
            total += sum(event.injection_rate) * dt
        end

        # If we've reached or passed 'time', stop
        if time <= t_end
            break
        end
    end
    return total
end

function get_all_parents(tstruct::TrapStructure, trap_id::Int)::Vector{Int}
    parents = Int[]
    current_id = trap_id
    while true
        parent_id = parentof(tstruct, current_id)

        if isnothing(parent_id)
            break
        end
        push!(parents, parent_id)
        current_id = parent_id
    end
    return parents
end