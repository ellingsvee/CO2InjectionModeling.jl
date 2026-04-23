import Interpolations
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, subtrapsof, FilledAmount
using Distributions: Normal, LogNormal, Uniform, truncated

export compute_leakage_volume, initialize_leakage_state, find_leakage_location
export compute_leakage_time_estimate, generate_leakage_weather_events
export get_true_topography_bottom
export compute_drainable_volume, compute_volume_with_drainage
export compute_dynamic_equilibrium_volume
export get_initial_leakage_time_estimates
"""
Stores the estimated time when a trap will reach its leakage volume. Similar to SWIM's ChangeTimeEstimate.
"""
struct LeakageTimeEstimate
    trap_id::Int
    min_time::Float64  # Earliest possible leakage time
    max_time::Float64  # Latest possible leakage time
end

"""
Find the grid cell where leakage occurs for a trap.
"""
function find_leakage_location(trap_id::Int, tstruct::TrapStructure)::CartesianIndex{2}
    footprint = tstruct.footprints[trap_id]
    topo_vals = tstruct.topography[footprint]
    min_idx = argmin(topo_vals)
    linear_idx = footprint[min_idx]
    return CartesianIndices(size(tstruct.topography))[linear_idx]
end


"""
Get the TRUE minimum topography elevation in a trap's footprint.
This is the actual lowest point in the footprint, including the area covered by child traps.
"""
function get_true_topography_bottom(trap_id::Int, tstruct::TrapStructure)::Float64
    footprint = tstruct.footprints[trap_id]
    return minimum(tstruct.topography[footprint])
end

"""
Compute the volume at which a trap reaches the leakage height threshold.

The CO2 column height is measured from the TRUE topography bottom (including child traps)
to the current water level. This correctly accounts for CO2 that fills child traps
before spilling to parent traps.

Returns:
- Float64: The volume at which CO2 height = leakage_height
- nothing: If the trap would spill before reaching leakage_height (leakage cannot occur)
"""
function compute_leakage_volume(
    trap_id::Int,
    z_vol_table::Tuple{Vector{Float64},Vector{Float64}},
    tstruct::TrapStructure,
    leakage_height::Float64
)::Union{Float64,Nothing}

    # Get the TRUE topography bottom (not the child spillpoint!)
    # This is the actual lowest point in the trap's footprint, including child areas
    true_bottom = get_true_topography_bottom(trap_id, tstruct)

    # Leakage occurs when water level reaches true_bottom + leakage_height
    leakage_elevation = true_bottom + leakage_height

    # Get spillpoint elevation (maximum fill level)
    spillpoint_elevation = tstruct.spillpoints[trap_id].elevation

    # If leakage elevation is above spillpoint, trap spills before it can leak
    if leakage_elevation >= spillpoint_elevation
        return nothing
    end

    # Use z_vol_table to find volume at leakage_elevation
    zvals, vvals = z_vol_table

    # Handle edge cases
    if length(zvals) == 1
        # Degenerate trap - no volume
        return 0.0
    end

    if leakage_elevation <= zvals[1]
        # Leakage elevation is below the trap's z_vol_table minimum.
        # For parent traps, zvals[1] is the child spillpoint.
        # This means leakage would occur while children are still filling,
        # so this trap's own volume at leakage is 0.
        return 0.0
    end

    if leakage_elevation >= zvals[end]
        # Leakage height is above trap capacity - shouldn't happen given check above
        return nothing
    end

    # Create interpolation function from z to volume
    z2v = Interpolations.linear_interpolation(zvals, vvals, extrapolation_bc=Interpolations.Line())

    return z2v(leakage_elevation)
end

"""
Compute the volume of CO2 that can drain from a trap during residual leakage.

The drainable volume is the portion that will leak out over the residual leakage time.
The residual (non-drainable) fraction remains trapped in pore spaces.

`residual_volume_fraction` should be `S_r / (1 - S_wi)` — the SWIM volume fraction
remaining after drainage, not the raw pore-space saturation `S_r`.
"""
function compute_drainable_volume(initial_volume::Float64, residual_volume_fraction::Float64)::Float64
    return initial_volume * (1.0 - residual_volume_fraction)
end


"""
Compute the current stored volume in a draining trap, accounting for residual drainage.
"""
function compute_volume_with_drainage(
    trap_id::Int,
    current_time::Float64,
    leakage_state::LeakageState
)::Union{Float64,Nothing}
    # If not draining, return nothing - caller should use actual volume from SWIM
    if !leakage_state.draining[trap_id]
        return nothing
    end

    # Get leakage parameters
    t_leak = leakage_state.leakage_start_time[trap_id]

    # If current_time is before leakage started, the trap wasn't leaking yet
    # Return nothing so caller uses actual volume from SWIM
    if current_time < t_leak
        return nothing
    end

    initial_vol = leakage_state.initial_volume_at_leak[trap_id]
    res_frac = leakage_state.residual_volume_fraction
    residual_time = leakage_state.residual_leakage_time

    # If no drainage (residual_time is Inf or 0, or res_frac is 1)
    if !isfinite(residual_time) || residual_time <= 0.0 || res_frac >= 1.0
        return initial_vol
    end

    # Compute time since leakage started
    time_since_leak = current_time - t_leak

    # Compute residual volume (what remains after full drainage)
    residual_vol = initial_vol * res_frac

    # Compute drainable volume
    drainable_vol = initial_vol - residual_vol

    # Linear drainage over residual_leakage_time
    if time_since_leak >= residual_time
        # Drainage complete
        return residual_vol
    else
        # Partial drainage
        fraction_drained = time_since_leak / residual_time
        return initial_vol - drainable_vol * fraction_drained
    end
end



"""
    compute_dynamic_equilibrium_volume(trap_id, current_time, leakage_state) -> Float64

Compute the stored volume in a trap accounting for leakage and drainage.

Directly leaking traps maintain capillary-gravity equilibrium: the CO2 column
stays at the threshold height (volume = leakage_volume) indefinitely. The
buoyancy pressure exactly balances the capillary entry pressure, so there is
no net driving force to push CO2 through the seal once inflow stops.

Non-leaking draining traps (descendants) experience residual drainage through
the spill hierarchy via `compute_volume_with_drainage`.
"""
function compute_dynamic_equilibrium_volume(
    trap_id::Int,
    current_time::Float64,
    leakage_state::LeakageState
)::Union{Float64,Nothing}

    if !leakage_state.leaking[trap_id]
        return compute_volume_with_drainage(trap_id, current_time, leakage_state)
    end

    t_leak = leakage_state.leakage_start_time[trap_id]
    if current_time < t_leak
        return nothing
    end

    return leakage_state.leakage_volume[trap_id]
end


"""
Initialize the leakage state for a layer, precomputing leakage volumes for all traps.

The `residual_volume_fraction` stored in LeakageState is computed as `S_r / (1 - S_wi)`,
converting from pore-space CO2 saturation to the equivalent fraction of SWIM volume.
This is because SWIM volumes represent pore space at maximum CO2 saturation `(1 - S_wi)`,
so the residual CO2 (at saturation `S_r`) corresponds to `S_r / (1 - S_wi)` of the SWIM volume.
"""
function initialize_leakage_state(
    tstruct::TrapStructure,
    z_vol_tables::Vector{Tuple{Vector{Float64},Vector{Float64}}},
    reservoir_properties::ReservoirProperties,
    residual_saturation::Float64,
    residual_leakage_time::Float64
)::LeakageState

    num_traps = numtraps(tstruct)

    if reservoir_properties.leakage_height isa Float64
        trap_leakage_heights = fill(reservoir_properties.leakage_height, num_traps)
    else
        trap_leakage_heights = reservoir_properties.leakage_height
    end

    # Precompute leakage volumes for all traps using their specific heights
    leakage_volumes = zeros(Float64, num_traps)
    for trap_id in 1:num_traps
        trap_height = trap_leakage_heights[trap_id]
        vol = compute_leakage_volume(trap_id, z_vol_tables[trap_id], tstruct, trap_height)
        # Use Inf for traps that cannot leak (spill before reaching leakage height)
        leakage_volumes[trap_id] = isnothing(vol) ? Inf : vol
    end

    # Correct parent traps that have shale-break children (leakage_height=0).
    # compute_leakage_volume measures column height from the TRUE topography
    # bottom of the entire footprint, which includes break children's deep
    # areas.  But a break child drains immediately, so the continuous CO2
    # column for the parent doesn't extend through that child's footprint.
    # Recompute using an effective bottom that excludes break children's areas.
    for trap_id in 1:num_traps
        trap_leakage_heights[trap_id] > 0.0 || continue

        children = subtrapsof(tstruct, trap_id)
        isempty(children) && continue

        break_children = [c for c in children if trap_leakage_heights[c] == 0.0]
        isempty(break_children) && continue

        # Collect all cells belonging to break children's footprints
        break_cells = Set{Int}()
        for bc in break_children
            union!(break_cells, tstruct.footprints[bc])
        end

        # Compute effective bottom from remaining cells
        parent_footprint = tstruct.footprints[trap_id]
        remaining_cells = [idx for idx in parent_footprint if !(idx in break_cells)]

        if isempty(remaining_cells)
            effective_bottom = z_vol_tables[trap_id][1][1]
        else
            effective_bottom = minimum(tstruct.topography[remaining_cells])
        end

        leakage_elevation = effective_bottom + trap_leakage_heights[trap_id]
        spillpoint_elevation = tstruct.spillpoints[trap_id].elevation

        if leakage_elevation >= spillpoint_elevation
            leakage_volumes[trap_id] = Inf
        else
            zvals, vvals = z_vol_tables[trap_id]
            if leakage_elevation <= zvals[1]
                leakage_volumes[trap_id] = 0.0
            elseif leakage_elevation >= zvals[end]
                leakage_volumes[trap_id] = Inf
            else
                z2v = Interpolations.linear_interpolation(
                    zvals, vvals, extrapolation_bc=Interpolations.Line())
                leakage_volumes[trap_id] = z2v(leakage_elevation)
            end
        end
    end

    # Convert raw S_r to the SWIM volume fraction: S_r / (1 - S_wi).
    # SWIM volumes are proportional to CO2 mass via porosity * (1 - S_wi) * cell_area.
    # After drainage, CO2 saturation drops from (1 - S_wi) to S_r, so the fraction
    # of SWIM volume (≡ CO2 mass) remaining is S_r / (1 - S_wi).
    S_wi = reservoir_properties.sand_irreducible_water_saturation
    residual_vol_fraction = residual_saturation / (1.0 - S_wi)

    return LeakageState(
        fill(false, num_traps),
        fill(false, num_traps),
        leakage_volumes,
        fill(Inf, num_traps),
        LeakageRecord[],
        trap_leakage_heights,
        fill(0.0, num_traps),
        residual_vol_fraction,
        residual_leakage_time,
    )
end




"""
Compute when a trap will reach its leakage threshold.
"""
function compute_leakage_time_estimate(
    trap_id::Int,
    cur_amounts::Vector{FilledAmount},
    cur_time::Float64,
    rateinfo,  # SWIM's RateInfo
    leakage_state::LeakageState,
    filled_traps::Vector{Bool},
    tstruct::TrapStructure
)::LeakageTimeEstimate

    # Already leaking - no future leakage event
    if leakage_state.leaking[trap_id]
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    # Get target volume for leakage
    target_vol = leakage_state.leakage_volume[trap_id]

    # If target is Inf, trap cannot leak (spills before reaching leakage height)
    if target_vol == Inf
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    # Get current volume
    current_vol = cur_amounts[trap_id].amount

    # If trap is not accumulating (children not filled), cannot estimate leakage
    children = subtrapsof(tstruct, trap_id)
    if !all(filled_traps[children])
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    # Check if already at or above target
    # Special case: if target_vol is 0, we need current_vol > 0 to trigger
    # (can't leak if there's no CO2 yet)
    if target_vol == 0.0
        if current_vol > 0.0
            # Trap has CO2 and should leak immediately (already above threshold)
            return LeakageTimeEstimate(trap_id, cur_time, cur_time)
        end
        # target_vol = 0 but no CO2 yet - will leak as soon as it gets any CO2
        # Continue to compute time estimate based on inflow
    elseif current_vol >= target_vol
        return LeakageTimeEstimate(trap_id, cur_time, cur_time)
    end

    # Compute time using inflow rate
    # Use SWIM's RateInfo accessors
    inflow = SurfaceWaterIntegratedModeling.getinflow(rateinfo, trap_id)
    smax = SurfaceWaterIntegratedModeling.getsmax(rateinfo, trap_id)
    smin = SurfaceWaterIntegratedModeling.getsmin(rateinfo, trap_id)

    # Net inflow bounds (same logic as SWIM's _compute_changetime_estimate)
    min_net_inflow = inflow - (smax - smin)
    max_net_inflow = inflow

    # If inflow is zero or negative, trap won't fill to leakage level
    if max_net_inflow <= 0
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    volume_needed = target_vol - current_vol

    # Time estimates (min time uses max inflow, max time uses min inflow)
    min_time = volume_needed / max_net_inflow
    max_time = (min_net_inflow > 0) ? volume_needed / min_net_inflow : Inf

    return LeakageTimeEstimate(
        trap_id,
        cur_time + min_time,
        cur_time + max_time
    )
end


"""
Initialize leakage time estimates for all traps.
"""
function get_initial_leakage_time_estimates(
    cur_amounts::Vector{FilledAmount},
    cur_time::Float64,
    rateinfo,
    leakage_state::LeakageState,
    filled_traps::Vector{Bool},
    tstruct::TrapStructure
)::Vector{LeakageTimeEstimate}

    return [
        compute_leakage_time_estimate(trap_id, cur_amounts, cur_time, rateinfo,
            leakage_state, filled_traps, tstruct)
        for trap_id in 1:numtraps(tstruct)
    ]
end


"""
Update leakage time estimates for affected traps.
"""
function update_leakage_time_estimates!(
    leakage_time_est::Vector{LeakageTimeEstimate},
    affected_traps::Vector{Int},
    cur_amounts::Vector{FilledAmount},
    cur_time::Float64,
    rateinfo,
    leakage_state::LeakageState,
    filled_traps::Vector{Bool},
    tstruct::TrapStructure
)
    for trap_id in affected_traps
        leakage_time_est[trap_id] = compute_leakage_time_estimate(
            trap_id, cur_amounts, cur_time, rateinfo,
            leakage_state, filled_traps, tstruct
        )
    end
end


"""
Find the next leakage event (minimum leakage time).
Returns (time, trap_id). If no leakage is pending, returns (Inf, 0).
"""
function find_next_leakage_event(leakage_time_est::Vector{LeakageTimeEstimate})::Tuple{Float64,Int}
    min_time = Inf
    min_trap = 0

    for est in leakage_time_est
        if est.min_time < min_time
            min_time = est.min_time
            min_trap = est.trap_id
        end
    end

    return (min_time, min_trap)
end


"""
Get the total injection rate at a given time from weather events.
"""
function _get_injection_rate_at_time(
    weather_events::Vector{WeatherEvent},
    t::Float64,
    grid_size::Tuple{Int,Int}
)::Float64
    # Find the weather event active at time t (last event with timestamp <= t)
    current_we = weather_events[1]
    for we in weather_events
        if we.timestamp <= t
            current_we = we
        else
            break
        end
    end
    if current_we.rain_rate isa Matrix
        return sum(current_we.rain_rate)
    else
        # Scalar rain_rate - uniform across domain
        return current_we.rain_rate * prod(grid_size)
    end
end

"""
Generate `WeatherEvent`s for an overlying layer that receives CO2 leaked from
the layer below.
"""
function generate_leakage_weather_events(
    seq::Vector{SpillEvent},
    leakage_state::LeakageState,
    tstruct::TrapStructure,
    rp_source::ReservoirProperties,
    rp_target::ReservoirProperties,
    direct_events::Vector{WeatherEvent};
    leakage_radius::Int=0,
    target_regions::Union{Nothing,Matrix{Int}}=nothing,
)::Vector{WeatherEvent}

    n_traps = numtraps(tstruct)
    grid_size = size(tstruct.topography)

    # Quick exit: nothing drains or leaks in the source layer.
    if !any(leakage_state.draining) && !any(leakage_state.leaking)
        return copy(direct_events)
    end

    # Unit conversion (includes density ratio for mass conservation across layers)
    density_ratio = rp_source.co2_density / rp_target.co2_density
    unit_factor = density_ratio * full_volume_to_rock_volume_scaling(rp_source) /
                  full_volume_to_rock_volume_scaling(rp_target)

    # Map each draining trap to its caprock exit location.
    drain_loc = Dict{Int,CartesianIndex{2}}()
    leak_locs = Dict{Int,CartesianIndex{2}}()
    for record in leakage_state.leakage_records
        trap = record.trap_id
        loc = record.leakage_location
        drain_loc[trap] = loc
        leak_locs[trap] = loc
        for desc in get_all_descendants(tstruct, trap)
            leakage_state.draining[desc] && (drain_loc[desc] = loc)
        end
    end

    # Also map draining traps that aren't descendants of a leaking trap.
    for trap in 1:n_traps
        leakage_state.draining[trap] || continue
        haskey(drain_loc, trap) && continue
        # Trace spill graph to find a trap with a known exit location
        t = tstruct.spillpoints[trap].downstream_region
        visited = Set{Int}(trap)
        while t > 0 && t <= n_traps && !(t in visited)
            if haskey(drain_loc, t)
                drain_loc[trap] = drain_loc[t]
                break
            end
            push!(visited, t)
            t = tstruct.spillpoints[t].downstream_region
        end
    end

    # Residual drainage parameters
    T_res = leakage_state.residual_leakage_time
    res_frac = leakage_state.residual_volume_fraction
    has_finite_drainage = isfinite(T_res) && T_res > 0.0 && res_frac < 1.0

    # Drainage rates for descendant (non-leaking) draining traps only.
    # Leaking traps maintain equilibrium and do not drain through the caprock.
    drain_rates = zeros(Float64, n_traps)
    if has_finite_drainage
        for trap in 1:n_traps
            leakage_state.leaking[trap] && continue
            leakage_state.draining[trap] || continue
            haskey(drain_loc, trap) || continue
            drain_rates[trap] =
                leakage_state.initial_volume_at_leak[trap] * (1 - res_frac) / T_res
        end
    end

    # Collect all timestamps where the combined rate can change
    timestamps = Set{Float64}()
    push!(timestamps, seq[1].timestamp)
    for se in seq
        push!(timestamps, se.timestamp)
    end
    for de in direct_events
        push!(timestamps, de.timestamp)
    end
    for trap in 1:n_traps
        t0 = leakage_state.leakage_start_time[trap]
        isfinite(t0) || continue
        push!(timestamps, t0)
        if has_finite_drainage && leakage_state.draining[trap] && !leakage_state.leaking[trap]
            push!(timestamps, t0 + T_res)
        end
    end

    sorted_ts = sort(collect(timestamps))

    # Build weather events
    result = WeatherEvent[]
    last_rain = nothing

    for (ti, t) in enumerate(sorted_ts)
        rain = zeros(Float64, grid_size)

        # Direct injection into the target layer
        if !isempty(direct_events)
            ix = findlast(de -> de.timestamp <= t, direct_events)
            if !isnothing(ix)
                rr = direct_events[ix].rain_rate
                rain .+= (rr isa Matrix ? rr : fill(rr, grid_size))
            end
        end

        seq_ix = findlast(se -> se.timestamp <= t, seq)

        # Leaking traps: all inflow passes through (IP equilibrium)
        for trap in 1:n_traps
            leakage_state.leaking[trap] || continue
            haskey(drain_loc, trap) || continue
            leakage_state.leakage_start_time[trap] > t && continue

            rate = isnothing(seq_ix) ? 0.0 : max(0.0, inflow_at(seq, seq_ix)[trap])

            if rate > 0
                spread_rate!(rain, Tuple(drain_loc[trap]), rate * unit_factor, leakage_radius; regions=target_regions)
            end
        end

        # Non-leaking draining traps: constant-rate drainage through spill hierarchy
        if has_finite_drainage
            for trap in 1:n_traps
                leakage_state.leaking[trap] && continue
                drain_rates[trap] == 0.0 && continue
                t0 = leakage_state.leakage_start_time[trap]
                (t >= t0 && t < t0 + T_res) || continue
                spread_rate!(rain, Tuple(drain_loc[trap]), drain_rates[trap] * unit_factor, leakage_radius; regions=target_regions)
            end
        end

        if isnothing(last_rain) || rain != last_rain
            push!(result, WeatherEvent(t, copy(rain)))
            last_rain = copy(rain)
        end
    end

    # Append zero-rate event at the latest drainage end time.
    if !isempty(result) && any(result[end].rain_rate .> 0.0)
        max_t_end = -Inf
        for trap in 1:n_traps
            t0 = leakage_state.leakage_start_time[trap]
            isfinite(t0) || continue
            if leakage_state.draining[trap] && !leakage_state.leaking[trap]
                max_t_end = max(max_t_end, t0 + T_res)
            end
        end
        if isfinite(max_t_end) && max_t_end > result[end].timestamp
            push!(result, WeatherEvent(max_t_end, zeros(Float64, grid_size)))
        end
    end

    return result
end
