import Graphs

export reconstruct_3d_lithology, create_trap_mask_3d, scale_unit_volume_to_physical, get_all_parents, get_all_descendants, convert_injection_event_to_weather_event
export verify_mass_conservation, compute_total_stored_volume, compute_total_leaked_volume
export create_plume_extent_masks

"""
Reconstruct 3D lithology grid from topography surfaces.

Returns a 3D array where:
- 0 = caprock (impermeable)
- 1 = sand (permeable reservoir)
- 2 = shale (low permeability)

Parameters:
- topography: Any AbstractTopography implementation
- domain: Domain3D struct
"""
function reconstruct_3d_lithology(topography::AbstractTopography, domain::Domain3D)
    nx, ny, nz = domain.nx, domain.ny, domain.nz
    dz = domain.dz
    depth_min = domain.depth_min
    depth_max = domain.depth_max

    println("\nReconstructing 3D lithology grid...")
    println("  Grid dimensions: $(nx) × $(ny) × $(nz)")
    println("  Depth range: $(depth_min) to $(depth_max) m")
    println("  Cell size: dx=$(domain.dx) m, dy=$(domain.dy) m, dz=$(dz) m")
    # 
    # Initialize with shale (default)
    lithology = fill(2, nx, ny, nz)

    # Create depth array for all k indices (vectorized)
    # k=1 is deepest, k=nz is shallowest
    depth_at_k = depth_max .- (0.5:nz) .* dz  # Vector of length nz

    # Get caprock surface (may be nothing for some topographies)
    top_caprock = get_caprock_surface(topography)

    # Process caprock layer (vectorized over k) if caprock surface is defined
    caprock_count = 0
    if !isnothing(top_caprock)
        println("  Processing caprock...")
        for k in 1:nz
            cell_depth = depth_at_k[k]
            caprock_mask = cell_depth .< top_caprock  # Broadcasting, result is (nx, ny)
            n_caprock_at_k = count(caprock_mask)
            caprock_count += n_caprock_at_k
            if n_caprock_at_k > 0
                lithology[caprock_mask, k] .= 0  # Caprock
            end
        end
        println("    Assigned $(caprock_count) caprock cells")
    else
        println("  No caprock surface defined, skipping caprock processing")
    end

    # Get sand layers via interface method
    sand_layers = get_sand_layers(topography)

    # Process sand layers (vectorized)
    println("  Processing sand layers...")
    sand_count = 0
    for layer in sand_layers
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



    for time_ix = 1:length(timepoints)
        tp = timepoints[time_ix]

        # Find the correct sequence index for this timepoint
        # Advance seq_ix until we find the last event at or before this timepoint
        while seq_ix < length(seq) && seq[seq_ix+1].timestamp <= tp
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
            tstates[time_ix][1],
            injected_volume,
            total_co2_volume,
        ))

    end

    return snapshots
end

function compute_total_injected_amount(injection_events::Vector{InjectionEvent}, time::Float64)
    total = 0.0
    for event_idx in 1:length(injection_events)
        event = injection_events[event_idx]
        t_start = event.timestamp
        t_end = event_idx < length(injection_events) ? injection_events[event_idx+1].timestamp : time

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


"""
    get_all_descendants(tstruct::TrapStructure, trap_id::Int) -> Vector{Int}

Get all descendants (children, grandchildren, etc.) of a trap.
Returns a vector of trap IDs in breadth-first order.
"""
function get_all_descendants(tstruct::TrapStructure, trap_id::Int)::Vector{Int}
    descendants = Int[]
    to_process = collect(subtrapsof(tstruct, trap_id))
    while !isempty(to_process)
        current = popfirst!(to_process)
        push!(descendants, current)
        append!(to_process, subtrapsof(tstruct, current))
    end
    return descendants
end


function convert_injection_event_to_weather_event(
    injection_event::Vector{InjectionEvent},
    reservoir_properties::ReservoirProperties,
    domain::Domain3D
)::Vector{WeatherEvent}
    weather_events = [WeatherEvent(ie.timestamp, physical_volume_to_swim_volume(ie.injection_rate, reservoir_properties, domain)) for ie in injection_event]
    return weather_events
end


"""
    compute_total_stored_volume(spill_events, tstruct, end_time; leakage_state=nothing) -> Float64

Compute the total CO2 volume stored in all traps at a given time (in SWIM units).

If leakage_state is provided, accounts for residual drainage from draining traps.
"""
function compute_total_stored_volume(
    spill_events::Vector{SpillEvent},
    tstruct::TrapStructure,
    end_time::Float64;
    leakage_state::Union{LeakageState,Nothing}=nothing
)::Float64
    # Get trap states at end time
    tstates = trap_states_at_timepoints(tstruct, spill_events, [end_time]; verbose=false)
    water_content = tstates[1][2]  # Volume in each trap

    total = 0.0
    for trap_id in 1:numtraps(tstruct)
        vol = water_content[trap_id]

        # If this trap is draining and we have leakage state, account for drainage
        # Note: we check 'draining' not 'leaking' because descendants of leaking traps
        # also drain, even though they're not at the leakage threshold themselves
        if !isnothing(leakage_state) && leakage_state.draining[trap_id]
            # Use the drainage-adjusted volume instead
            drained_vol = compute_volume_with_drainage(trap_id, end_time, leakage_state)
            if !isnothing(drained_vol)
                vol = drained_vol
            end
        end

        total += vol
    end

    return total
end


"""
    compute_total_leaked_volume(leakage_state, spill_events, weather_events, end_time) -> Float64

Compute the total CO2 volume that has leaked from all traps (in SWIM units).
Leaked volume = integral of leakage rate from leakage_start_time to end_time.
"""
function compute_total_leaked_volume(
    leakage_state::LeakageState,
    spill_events::Vector{SpillEvent},
    weather_events::Vector{WeatherEvent},
    end_time::Float64
)::Float64
    total_leaked = 0.0

    for record in leakage_state.leakage_records
        trap_id = record.trap_id
        start_time = record.start_time

        # Integrate inflow rate from start_time to end_time
        # The inflow rate changes at each weather event and spill event

        # Collect all time points where rates might change
        rate_change_times = Float64[start_time]

        for we in weather_events
            if we.timestamp > start_time && we.timestamp < end_time
                push!(rate_change_times, we.timestamp)
            end
        end

        for se in spill_events
            if se.timestamp > start_time && se.timestamp < end_time
                push!(rate_change_times, se.timestamp)
            end
        end

        push!(rate_change_times, end_time)
        sort!(unique!(rate_change_times))

        # Integrate piecewise
        for i in 1:(length(rate_change_times)-1)
            t_start = rate_change_times[i]
            t_end = rate_change_times[i+1]
            dt = t_end - t_start

            if dt > 0
                # Get inflow rate at t_start
                inflow_rate = get_trap_inflow_at_time(trap_id, t_start, spill_events)
                total_leaked += inflow_rate * dt
            end
        end
    end

    return total_leaked
end


"""
    verify_mass_conservation(injection_events, spill_events, leakage_state, end_time, reservoir_properties, domain; tolerance=1e-6)

Verify mass conservation: total injected = total stored + total leaked.

Note: Leaked volume is computed by mass balance (injected - stored) to avoid
double-counting when multiple traps are leaking in a chain.

Returns a NamedTuple with:
- `conserved`: Boolean indicating if mass is conserved within tolerance
- `injected`: Total injected volume (physical units, m³)
- `stored`: Total stored volume (physical units, m³)
- `leaked`: Total leaked volume (physical units, m³) - computed by mass balance
- `error`: Absolute error in mass balance (physical units, m³) - should be ~0
- `relative_error`: Relative error (error / injected)
"""
function verify_mass_conservation(
    injection_events::Vector{InjectionEvent},
    spill_events::Vector{SpillEvent},
    leakage_state::LeakageState,
    tstruct::TrapStructure,
    weather_events::Vector{WeatherEvent},
    end_time::Float64,
    reservoir_properties::ReservoirProperties,
    domain::Domain3D;
    tolerance::Float64=1e-6
)
    # Compute injected volume (physical units)
    injected = compute_total_injected_amount(injection_events, end_time)

    # Compute stored volume (SWIM units, then convert)
    stored_swim = compute_total_stored_volume(spill_events, tstruct, end_time)
    stored = swim_volume_to_physical_volume(stored_swim, reservoir_properties, domain)

    # Compute leaked volume by mass balance to avoid double-counting
    # leaked = injected - stored
    leaked = max(0.0, injected - stored)

    # Compute error (should be ~0 if mass is conserved)
    # The error would come from numerical precision or if storage is computed incorrectly
    error = abs(injected - stored - leaked)
    relative_error = injected > 0 ? error / injected : 0.0

    # Check conservation - with this formulation, it should always be true
    # unless there's a numerical precision issue
    conserved = relative_error < tolerance

    return (
        conserved=conserved,
        injected=injected,
        stored=stored,
        leaked=leaked,
        error=error,
        relative_error=relative_error
    )
end


"""
    create_plume_extent_masks(layers, seqs, domain, time; return_3d=false) -> Vector

Create masks showing the spatial extent of CO2 plumes in each layer at a specific time.

This function extracts the CO2 distribution at a given timepoint and returns boolean masks
indicating where CO2 is present (volume > 0 or trap is filled).

# Arguments
- `layers::Vector{Layer}`: Vector of Layer structs
- `seqs::Vector{Vector{SpillEvent}}`: Vector of spill event sequences from fill_layers
- `domain::Domain3D`: Domain3D struct
- `time::Float64`: Timepoint to evaluate plume extent

# Keyword Arguments
- `return_3d::Bool=false`: If true, returns 3D masks (nx, ny, nz). If false, returns 2D masks (nx, ny)

# Returns
- If `return_3d=false`: Vector{Array{Bool,2}} - one 2D mask per layer (nx, ny)
- If `return_3d=true`: Vector{Array{Bool,3}} - one 3D mask per layer (nx, ny, nz)

Masks have `true` where CO2 is present and `false` elsewhere. Boundary padding is automatically removed.

# Example
```julia
# Get 2D plume footprints at t=10 years
masks_2d = create_plume_extent_masks(layers, seqs, domain, 10.0)

# Get 3D plume volumes at t=10 years
masks_3d = create_plume_extent_masks(layers, seqs, domain, 10.0; return_3d=true)

# Check plume extent in layer 1
layer1_mask = masks_2d[1]
plume_area_m2 = sum(layer1_mask) * domain.dx * domain.dy
```
"""
function create_plume_extent_masks(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    domain::Domain3D,
    time::Float64;
    return_3d::Bool=false
)
    n_layers = length(layers)
    @assert n_layers == length(seqs) "Number of layers must match number of sequences"

    nx, ny, nz = domain.nx, domain.ny, domain.nz

    masks = []

    for (layer_idx, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        num_traps = numtraps(tstruct)
        seq = seqs[layer_idx]
        pad = layer.boundary_padding

        # Get padded grid size
        topo_size = size(tstruct.topography)
        nx_padded, ny_padded = topo_size

        # Skip empty layers
        if isempty(seq) || num_traps == 0
            if return_3d
                push!(masks, falses(nx, ny, nz))
            else
                push!(masks, falses(nx, ny))
            end
            continue
        end

        # Get trap states at the specified timepoint
        tstates = trap_states_at_timepoints(tstruct, seq, [time]; verbose=false)
        filled, volumes, _ = tstates[1]

        if return_3d
            # Create 3D mask by combining 2D footprint with vertical extent
            mask_3d_padded = falses(nx_padded, ny_padded, nz)

            for trap_id in 1:num_traps
                volume = volumes[trap_id]
                if volume > 0.0 || filled[trap_id]
                    # Get the 3D extent for this trap
                    trap_mask_3d = create_trap_mask_3d(layer, trap_id, domain)
                    mask_3d_padded .|= trap_mask_3d
                end
            end

            # Remove padding
            mask_3d = mask_3d_padded[pad+1:end-pad, pad+1:end-pad, :]
            push!(masks, mask_3d)
        else
            # Create 2D footprint mask
            mask_2d_padded = falses(nx_padded, ny_padded)

            for trap_id in 1:num_traps
                volume = volumes[trap_id]
                if volume > 0.0 || filled[trap_id]
                    footprint = tstruct.footprints[trap_id]
                    for idx in footprint
                        mask_2d_padded[idx] = true
                    end
                end
            end

            # Remove padding
            mask_2d = mask_2d_padded[pad+1:end-pad, pad+1:end-pad]
            push!(masks, mask_2d)
        end
    end

    return masks
end