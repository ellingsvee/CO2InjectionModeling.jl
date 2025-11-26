export LeakageState, LeakageEvent
export compute_leakage_height, _compute_trap_height, compute_per_trap_leaked_volumes,  get_total_leaked_volume


struct LeakageEvent
    timestamp::Float64           # Time when leakage was detected (at spill event)
    trap_id::Int                # Trap that leaked
    height_at_detection::Float64  # Height when leakage was detected (>= threshold)
    volume_in_trap::Float64     # Volume in trap at detection (swim units)
    source_regions::Set{Int}    # Regions whose injection contributed to this leakage
    connected_filled_traps::Vector{Int}  # Trap IDs that are filled and connected (will leak)
    total_leakable_volume::Float64       # Total volume that will leak (after residual trapping)
    residual_trapped_volume::Float64     # Volume that stays trapped (residual)
end

mutable struct LeakageState
    leaked_traps::Vector{Bool}              # Which traps have leaked
    leakage_events::Vector{LeakageEvent}    # Detailed leakage events
    leakage_threshold::Float64              # Height threshold (meters)
    first_leakage_time::Union{Float64, Nothing}  # Time of first leakage (Nothing if no leakage yet)
end

"""
Helper function to compute the height of CO2 in a specific trap at the current state.
Returns the maximum height within the trap's footprint.
"""
function _compute_trap_height(trap_id::Int, cur_amounts::Vector, z_vol_tables, tstruct)
    water_volume = cur_amounts[trap_id].amount

    # If trap is empty, height is zero
    if water_volume <= 0.0
        return 0.0
    end

    # Get water elevation from volume
    water_elevation = _volume_to_elevation(water_volume, z_vol_tables[trap_id])

    min_base_elevation = minimum(tstruct.topography[tstruct.footprints[trap_id]])
    distance_from_spillpoint_to_trap_bottom = 0.0
    children = subtrapsof(tstruct, trap_id)
    if !isempty(children)
        distance_from_spillpoint_to_trap_bottom = max.(distance_from_spillpoint_to_trap_bottom, tstruct.spillpoints[children[1]].elevation - min_base_elevation)
    end

    water_elevation += distance_from_spillpoint_to_trap_bottom

    # Height is elevation difference
    return max(0.0, water_elevation - min_base_elevation)
end

"""
Compute the leakage height based on the reservoir properties.
Uses the relation: cappilary_pressure_threshold = density_difference_between_brine_and_co2 * g * column_height
"""
function compute_leakage_height(reservoir_properties::ReservoirProperties) 
    g = 9.81 # m/s^2
    height = reservoir_properties.shale_pressure_threshold / (reservoir_properties.brine_co2_density_difference * g)
    # height = 15.0 # TEMPORARY FIX FOR TESTING. Know leakage should occur here.
    return height
end

"""
Apply leakage state to rain_rate, zeroing out injection in regions that have already leaked.

This function is called at the start of each weather event to ensure that regions
which had their injection stopped due to leakage remain stopped across weather event boundaries.
"""
function _apply_leakage_to_rain_rate!(rain_rate::Matrix{Float64}, leakage_state::LeakageState,
                                      tstruct::TrapStructure)
    # For each leakage event, turn off injection in the contributing source regions
    for event in leakage_state.leakage_events
        for source_region in event.source_regions
            region_mask = tstruct.regions .== source_region
            rain_rate[region_mask] .= 0.0
        end
    end
end

"""
Trigger leakage for a trap and propagate to all parent traps.
Modifies leakage_state and rain_rate in place.

The leaked volume is only recorded for the trap that triggered the leakage (trap_id).
Parent traps are marked as leaked to stop injection, but their volumes are not added
to the total leaked amount (to avoid double counting).

Uses source_tracker to determine which injection sources contributed CO2 to this trap,
and turns off injection only in those source regions.

Computes the total volume that will leak based on:
- Only filled parent traps connected to the leaking trap contribute
- Only (1 - sand_residual_co2_saturation) fraction of the volume leaks
- The remaining fraction stays as residually trapped CO2
"""
function _trigger_leakage!(trap_id::Int, leakage_state::LeakageState,
                          cur_amounts::Vector, rain_rate::Matrix{Float64},
                          z_vol_tables, tstruct, cur_time::Float64,
                          source_tracker::SourceTracker, filled_traps::Vector{Bool},
                          reservoir_properties::ReservoirProperties, verbose::Bool)
    # Get all traps that need to be marked as leaked (trap + all parents)
    all_affected_traps = _identify_all_parent_traps(trap_id, tstruct)

    # Record first leakage time if this is the first leakage
    if isnothing(leakage_state.first_leakage_time)
        leakage_state.first_leakage_time = cur_time
    end

    # Find which injection sources contributed to this leakage
    contributing_sources = find_regions_to_stop(source_tracker, trap_id)

    # First, record the leakage event for the trap that triggered it
    if !leakage_state.leaked_traps[trap_id]
        volume_in_trap = cur_amounts[trap_id].amount
        # Since cur_time is the exact time when the threshold was reached,
        # the height at that time should be exactly the threshold
        height_at_detection = leakage_state.leakage_threshold

        # Find connected filled traps (these will contribute to leakage)
        connected_filled_traps = _find_connected_filled_parents(trap_id, filled_traps, tstruct)

        # Calculate total volume in connected filled traps
        total_connected_volume = 0.0
        for connected_trap in connected_filled_traps
            total_connected_volume += cur_amounts[connected_trap].amount
        end

        # Calculate leakable vs residual volumes
        # Only (1 - sand_residual_co2_saturation) fraction leaks through
        leakage_fraction = 1.0 - reservoir_properties.sand_residual_co2_saturation
        total_leakable_volume = total_connected_volume * leakage_fraction
        residual_trapped_volume = total_connected_volume * reservoir_properties.sand_residual_co2_saturation

        leakage_state.leaked_traps[trap_id] = true

        push!(leakage_state.leakage_events,
              LeakageEvent(cur_time, trap_id, height_at_detection, volume_in_trap,
                          copy(contributing_sources), connected_filled_traps,
                          total_leakable_volume, residual_trapped_volume))

        if verbose
            println("  Connected filled traps: $connected_filled_traps")
            println("  Total connected volume: $total_connected_volume")
            println("  Leakable volume: $total_leakable_volume ($(leakage_fraction*100)%)")
            println("  Residual trapped: $residual_trapped_volume ($(reservoir_properties.sand_residual_co2_saturation*100)%)")
        end
    end

    # Now mark all affected traps as leaked
    for affected_trap in all_affected_traps
        # Skip if already leaked
        if leakage_state.leaked_traps[affected_trap]
            continue
        end

        # Mark trap as leaked (but don't count volume again)
        leakage_state.leaked_traps[affected_trap] = true
    end

    if verbose
        println("  Turning off injection in source regions: $contributing_sources")
    end

    # Zero out injection in all regions that contributed to this leakage
    for source_region in contributing_sources
        region_mask = tstruct.regions .== source_region
        rain_rate[region_mask] .= 0.0
    end
end


function find_leakage_location(trap_id::Int, tstruct::TrapStructure)
    # Get the exact leakage location
    # This is the location in the trap footprint with the smallest topography elevation
    footprint = tstruct.footprints[trap_id]                # Vector of linear indices
    topography = tstruct.topography[footprint]             # Elevations at those indices
    min_idx = argmin(topography)                           # Index of minimum elevation within the footprint
    leakage_linear_idx = footprint[min_idx]                # Linear index in the domain
    leakage_cartesian_idx = CartesianIndices(size(tstruct.topography))[leakage_linear_idx]  # (i, j) coordinates
    return leakage_cartesian_idx
end


"""
Helper function to compute the exact time when a trap reaches a specific height threshold.

This function uses the trap's inflow rate and volume-height relationship to determine
precisely when the height threshold is crossed, rather than detecting it at discrete
trap fill events.

Returns:
- exact_time: The time when the threshold is reached (or nothing if not reached)
"""
function _compute_time_to_height_threshold(trap_id::Int, threshold_height::Float64,
                                           cur_amounts::Vector, z_vol_tables, tstruct,
                                           rateinfo, cur_time::Float64)
    # Get current volume and time
    current_volume = cur_amounts[trap_id].amount
    last_update_time = cur_amounts[trap_id].time

    # Compute current height
    current_height = _compute_trap_height(trap_id, cur_amounts, z_vol_tables, tstruct)

    # If already above threshold, return current time
    if current_height >= threshold_height
        return cur_time
    end

    # Get the minimum base elevation of the trap
    min_base_elevation = minimum(tstruct.topography[tstruct.footprints[trap_id]])

    # Adjust for parent trap geometry if needed
    distance_from_spillpoint_to_trap_bottom = 0.0
    children = SurfaceWaterIntegratedModeling.subtrapsof(tstruct, trap_id)
    if !isempty(children)
        distance_from_spillpoint_to_trap_bottom = max(distance_from_spillpoint_to_trap_bottom,
                                                       tstruct.spillpoints[children[1]].elevation - min_base_elevation)
    end

    # Compute target elevation from target height
    target_elevation = min_base_elevation + threshold_height - distance_from_spillpoint_to_trap_bottom

    # Convert target elevation to target volume using inverse interpolation
    z_vol_table = z_vol_tables[trap_id]
    zvals = z_vol_table[1]
    vvals = z_vol_table[2]

    # Check if target elevation is reachable
    if target_elevation < minimum(zvals) || target_elevation > maximum(zvals)
        return nothing  # Height not reachable in this trap
    end

    # Interpolate to find target volume
    target_volume = Interpolations.linear_interpolation(zvals, vvals, extrapolation_bc=Interpolations.Flat())(target_elevation)

    # If target volume is less than current, threshold was already passed
    if target_volume <= current_volume
        return cur_time
    end

    # Get inflow rate for this trap
    inflow_rate = SurfaceWaterIntegratedModeling.getinflow(rateinfo, trap_id)

    # If inflow is zero or negative, height won't be reached
    if inflow_rate <= 0.0
        return nothing
    end

    # Compute time needed to reach target volume from last update
    volume_needed = target_volume - current_volume
    time_from_last_update = volume_needed / inflow_rate

    # Return exact time when threshold is reached
    return last_update_time + time_from_last_update
end

"""
Check all traps for leakage condition and trigger leakage if needed.
Returns true if any leakage was detected and triggered.
"""
function _check_and_trigger_leakage!(leakage_state::LeakageState, cur_amounts::Vector,
                                     rain_rate::Matrix{Float64}, z_vol_tables, tstruct,
                                     cur_time::Float64, source_tracker::SourceTracker,
                                     filled_traps::Vector{Bool}, reservoir_properties::ReservoirProperties,
                                     verbose::Bool, rate_info=nothing)
    leakage_occurred = false
    num_traps = length(cur_amounts)

    for trap_id in 1:num_traps
        # Skip if already leaked
        if leakage_state.leaked_traps[trap_id]
            continue
        end

        # Compute current height
        height = _compute_trap_height(trap_id, cur_amounts, z_vol_tables, tstruct)

        # Check if exceeds threshold
        if height > leakage_state.leakage_threshold
            # Compute exact time when threshold was reached
            exact_leakage_time = cur_time
            if !isnothing(rate_info)
                computed_time = _compute_time_to_height_threshold(
                    trap_id, leakage_state.leakage_threshold,
                    cur_amounts, z_vol_tables, tstruct, rate_info, cur_time)
                if !isnothing(computed_time)
                    exact_leakage_time = computed_time
                end
            end

            if verbose
                println("Leakage detected in trap $trap_id at time $exact_leakage_time (height: $height m at detection, threshold: $(leakage_state.leakage_threshold) m)")
            end

            _trigger_leakage!(trap_id, leakage_state, cur_amounts, rain_rate,
                            z_vol_tables, tstruct, exact_leakage_time, source_tracker,
                            filled_traps, reservoir_properties, verbose)
            leakage_occurred = true
        end
    end

    return leakage_occurred
end

"""
Get total leaked volume (sum across all traps).

Parameters:
- leaked_volumes_per_trap: Dictionary from compute_per_trap_leaked_volumes()

Returns:
- Total leaked volume in swim units
"""
function get_total_leaked_volume(leaked_volumes_per_trap::Dict{Int, Float64})
    return sum(values(leaked_volumes_per_trap))
end