import Interpolations
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, subtrapsof, FilledAmount
using Distributions: Normal, LogNormal, Uniform, truncated

export compute_leakage_volume, initialize_leakage_state, find_leakage_location
export compute_leakage_time_estimate, generate_leakage_weather_events
export get_true_topography_bottom
export compute_drainable_volume, compute_volume_with_drainage, compute_residual_drainage_rate
export compute_dynamic_equilibrium_volume, update_leaking_trap_inflow_state!
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
"""
function compute_drainable_volume(initial_volume::Float64, residual_saturation::Float64)::Float64
    return initial_volume * (1.0 - residual_saturation)
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
    residual_sat = leakage_state.residual_saturation
    residual_time = leakage_state.residual_leakage_time

    # If no drainage (residual_time is Inf or 0, or residual_sat is 1)
    if !isfinite(residual_time) || residual_time <= 0.0 || residual_sat >= 1.0
        return initial_vol
    end

    # Compute time since leakage started
    time_since_leak = current_time - t_leak

    # Compute residual volume (what remains after full drainage)
    residual_vol = initial_vol * residual_sat

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
Compute the current residual drainage rate from a leaking trap.

This is the rate at which CO2 is draining from the trap due to residual leakage
(separate from the pass-through of new injections).
"""
function compute_residual_drainage_rate(
    trap_id::Int,
    current_time::Float64,
    leakage_state::LeakageState
)::Float64
    # If not leaking, no drainage
    if !leakage_state.leaking[trap_id]
        return 0.0
    end

    # Get leakage parameters
    t_leak = leakage_state.leakage_start_time[trap_id]
    initial_vol = leakage_state.initial_volume_at_leak[trap_id]
    residual_sat = leakage_state.residual_saturation
    residual_time = leakage_state.residual_leakage_time

    # If no drainage configured
    if !isfinite(residual_time) || residual_time <= 0.0 || residual_sat >= 1.0
        return 0.0
    end

    # Check if we're still in the drainage period
    time_since_leak = current_time - t_leak
    if time_since_leak < 0.0 || time_since_leak >= residual_time
        # Not yet started or already completed
        return 0.0
    end

    # Constant drainage rate during the drainage period
    drainable_vol = initial_vol * (1.0 - residual_sat)
    return drainable_vol / residual_time
end


"""
    compute_dynamic_equilibrium_volume(trap_id, current_time, leakage_state) -> Float64

Compute the stored volume in a directly leaking trap, accounting for dynamic
equilibrium. This means that rainage only accumulates during periods without 
inflow.

Due to the equilibrium assumption of IP theory, a leaking traps CO2 column is 
maintained at the threshold height (volume = leakage_volume) when it has positive inflow. 
Drainage only occurs when inflow drops to zero, representing buoyancy-driven flow through 
the caprock once the column is no longer sustained by injection.

For non-leaking draining traps (descendants), falls back to the standard
`compute_volume_with_drainage`.
"""
function compute_dynamic_equilibrium_volume(
    trap_id::Int,
    current_time::Float64,
    leakage_state::LeakageState
)::Union{Float64,Nothing}

    # For non-leaking draining traps, use standard drainage
    if !leakage_state.leaking[trap_id]
        return compute_volume_with_drainage(trap_id, current_time, leakage_state)
    end

    # Not yet leaking
    if !leakage_state.draining[trap_id] && !leakage_state.leaking[trap_id]
        return nothing
    end

    t_leak = leakage_state.leakage_start_time[trap_id]
    if current_time < t_leak
        return nothing
    end

    residual_sat = leakage_state.residual_saturation
    residual_time = leakage_state.residual_leakage_time

    V_leak = leakage_state.leakage_volume[trap_id]

    # If no drainage configured, volume stays at equilibrium
    if !isfinite(residual_time) || residual_time <= 0.0 || residual_sat >= 1.0
        return V_leak
    end

    # Use V_leak (not initial_volume_at_leak) as the drainage reference for
    # leaking traps. The equilibrium volume is V_leak; any overshoot beyond
    # V_leak at the moment leakage was detected is a discrete-timestep artifact.
    V_res = V_leak * residual_sat
    drainable_vol = V_leak * (1.0 - residual_sat)
    drain_rate = drainable_vol / residual_time

    t_last_change = leakage_state.time_of_last_state_change[trap_id]

    if current_time < t_last_change
        # Query time is before the last state transition.
        # The trap was in the OPPOSITE inflow state at this time.
        if !leakage_state.has_inflow[trap_id]
            # Currently no inflow, but at query time the trap HAD inflow.
            # Column was maintained at V_leak.
            return V_leak
        else
            # Currently has inflow, but at query time the trap had NO inflow (was draining).
            # Reconstruct volume: at the transition (t_last_change), volume was
            # volume_at_last_state_change. Going backward in time, volume increases
            # (un-draining) at drain_rate.
            V_at_transition = leakage_state.volume_at_last_state_change[trap_id]
            return min(V_leak, max(V_res, V_at_transition + drain_rate * (t_last_change - current_time)))
        end
    end

    if leakage_state.has_inflow[trap_id]
        # Trap has inflow — column maintained at V_leak
        return V_leak
    else
        # No inflow — draining from volume at last state change
        V_at_change = leakage_state.volume_at_last_state_change[trap_id]
        t_since_change = current_time - t_last_change
        drained = drain_rate * t_since_change
        return max(V_res, V_at_change - drained)
    end
end


"""
    update_leaking_trap_inflow_state!(leakage_state, trap_id, new_has_inflow, current_time)

Update the dynamic equilibrium tracking for a leaking trap when its inflow state changes.
Called from fill_layer.jl after flow rates are recomputed.
"""
function update_leaking_trap_inflow_state!(
    leakage_state::LeakageState,
    trap_id::Int,
    new_has_inflow::Bool,
    current_time::Float64
)
    !leakage_state.leaking[trap_id] && return

    old_has_inflow = leakage_state.has_inflow[trap_id]
    t_last = leakage_state.time_of_last_state_change[trap_id]

    if old_has_inflow == new_has_inflow
        return  # No state change
    end

    residual_sat = leakage_state.residual_saturation
    residual_time = leakage_state.residual_leakage_time
    has_finite_drainage = isfinite(residual_time) && residual_time > 0.0 && residual_sat < 1.0

    V_leak = leakage_state.leakage_volume[trap_id]
    V_res = V_leak * residual_sat

    if old_has_inflow && !new_has_inflow
        # Transition: inflow → no inflow. Column was maintained, now drainage begins.
        leakage_state.volume_at_last_state_change[trap_id] = V_leak
        leakage_state.time_of_last_state_change[trap_id] = current_time
        leakage_state.has_inflow[trap_id] = false

    elseif !old_has_inflow && new_has_inflow
        # Transition: no inflow → inflow. Compute volume after drainage, begin refilling.
        if has_finite_drainage && isfinite(t_last)
            dt = current_time - t_last
            drain_rate = V_leak * (1.0 - residual_sat) / residual_time
            V_start = leakage_state.volume_at_last_state_change[trap_id]
            V_current = max(V_res, V_start - drain_rate * dt)
            leakage_state.cumulative_no_inflow_time[trap_id] += dt
            leakage_state.volume_at_last_state_change[trap_id] = V_current
        else
            leakage_state.volume_at_last_state_change[trap_id] = V_leak
        end
        leakage_state.time_of_last_state_change[trap_id] = current_time
        leakage_state.has_inflow[trap_id] = true
    end
end


"""
Initialize the leakage state for a layer, precomputing leakage volumes for all traps.
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

    return LeakageState(
        fill(false, num_traps),           # leaking
        fill(false, num_traps),           # draining
        leakage_volumes,                   # leakage_volume
        fill(Inf, num_traps),             # leakage_start_time
        LeakageRecord[],                   # leakage_records
        trap_leakage_heights,              # leakage_height (NOW PER-TRAP VECTOR)
        fill(0.0, num_traps),             # initial_volume_at_leak (0 until leakage starts)
        residual_saturation,               # residual_saturation
        residual_leakage_time,             # residual_leakage_time
        fill(0.0, num_traps),             # cumulative_no_inflow_time
        fill(0.0, num_traps),             # volume_at_last_state_change
        fill(Inf, num_traps),             # time_of_last_state_change
        fill(false, num_traps)            # has_inflow
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
    sat = leakage_state.residual_saturation
    has_finite_drainage = isfinite(T_res) && T_res > 0.0 && sat < 1.0

    # Drainage rate for each draining trap (source-layer SWIM units).
    # For directly leaking traps: use V_leak as reference (not initial_volume_at_leak)
    #   to be consistent with compute_dynamic_equilibrium_volume and compute_total_drained.
    # For descendant (non-leaking) draining traps: use initial_volume_at_leak.
    drain_rates = zeros(Float64, n_traps)
    if has_finite_drainage
        for trap in 1:n_traps
            leakage_state.draining[trap] || continue
            haskey(drain_loc, trap) || continue
            if leakage_state.leaking[trap]
                drain_rates[trap] =
                    leakage_state.leakage_volume[trap] * (1 - sat) / T_res
            else
                drain_rates[trap] =
                    leakage_state.initial_volume_at_leak[trap] * (1 - sat) / T_res
            end
        end
    end

    # Dynamic equilibrium state for directly leaking traps
    # Track stored volume to determine whether trap is in equilibrium,
    # draining, or refilling at each timestamp.
    # Start at V_leak (equilibrium volume), not initial_volume_at_leak.
    trap_volume = zeros(Float64, n_traps)
    for trap in 1:n_traps
        leakage_state.leaking[trap] || continue
        trap_volume[trap] = leakage_state.leakage_volume[trap]
    end

    # Collect all timestamps where the combined rate can change
    timestamps = Set{Float64}()
    push!(timestamps, seq[1].timestamp)          # always start at sim beginning

    for se in seq
        push!(timestamps, se.timestamp)
    end  # passthrough rate changes
    for de in direct_events
        push!(timestamps, de.timestamp)
    end  # direct injection changes

    for trap in 1:n_traps
        t0 = leakage_state.leakage_start_time[trap]
        isfinite(t0) || continue
        push!(timestamps, t0)                        # leakage begins
        # Only add drainage end time for non-leaking draining traps (descendants).
        # Leaking traps use dynamic equilibrium, so their drainage timing depends
        # on inflow history and is computed below.
        if has_finite_drainage && leakage_state.draining[trap] && !leakage_state.leaking[trap]
            push!(timestamps, t0 + T_res)
        end
    end

    sorted_ts = sort(collect(timestamps))

    # Get inflow to a leaking trap at spill sequence index
    function _get_trap_inflow(trap, seq_ix)
        isnothing(seq_ix) && return 0.0
        inflows = inflow_at(seq, seq_ix)
        return max(0.0, inflows[trap])
    end

    # Compute dynamic equilibrium timeline for leaking traps
    extra_timestamps = Float64[]
    prev_inflow = Dict{Int,Float64}()  # previous inflow per leaking trap

    for (ti, t) in enumerate(sorted_ts)
        seq_ix = findlast(se -> se.timestamp <= t, seq)
        t_next = ti < length(sorted_ts) ? sorted_ts[ti+1] : Inf

        for trap in 1:n_traps
            leakage_state.leaking[trap] || continue
            haskey(drain_loc, trap) || continue
            leakage_state.leakage_start_time[trap] > t && continue

            V_leak = leakage_state.leakage_volume[trap]
            V_res = V_leak * sat
            dr = drain_rates[trap]

            inflow = _get_trap_inflow(trap, seq_ix)

            # Advance volume from previous interval.
            if ti > 1 && sorted_ts[ti-1] >= leakage_state.leakage_start_time[trap]
                dt = t - sorted_ts[ti-1]
                prev_q = get(prev_inflow, trap, 0.0)
                if dt > 0
                    if prev_q > 0 && trap_volume[trap] < V_leak
                        # Was refilling
                        trap_volume[trap] = min(trap_volume[trap] + prev_q * dt, V_leak)
                    elseif prev_q <= 0 && trap_volume[trap] > V_res && has_finite_drainage
                        # Was draining
                        trap_volume[trap] = max(trap_volume[trap] - dr * dt, V_res)
                    end
                    # Equilibrium (prev_q > 0 && V >= V_leak): no volume change
                end
            end

            prev_inflow[trap] = inflow

            # Check for mid-interval transitions in [t, t_next)
            if isfinite(t_next) && t_next > t
                dt_interval = t_next - t
                if inflow > 0 && trap_volume[trap] < V_leak
                    # Refilling: when does V reach V_leak?
                    time_to_refill = (V_leak - trap_volume[trap]) / inflow
                    if time_to_refill < dt_interval - 1e-12
                        push!(extra_timestamps, t + time_to_refill)
                    end
                elseif inflow <= 0 && trap_volume[trap] > V_res && has_finite_drainage && dr > 0
                    # Draining: when does V reach V_res?
                    time_to_drain = (trap_volume[trap] - V_res) / dr
                    if time_to_drain < dt_interval - 1e-12
                        push!(extra_timestamps, t + time_to_drain)
                    end
                end
            end
        end
    end

    # If we found mid-interval transitions, merge extra timestamps
    if !isempty(extra_timestamps)
        for et in extra_timestamps
            push!(timestamps, et)
        end
        sorted_ts = sort(collect(timestamps))
    end

    # Always reset trap volumes and prev_inflow before the main pass,
    # since the first pass consumes them.
    for trap in 1:n_traps
        leakage_state.leaking[trap] || continue
        trap_volume[trap] = leakage_state.leakage_volume[trap]
    end
    empty!(prev_inflow)

    # Build weather events
    result = WeatherEvent[]
    last_rain = nothing

    for (ti, t) in enumerate(sorted_ts)
        rain = zeros(Float64, grid_size)

        # Direct injection into the target layer.
        if !isempty(direct_events)
            ix = findlast(de -> de.timestamp <= t, direct_events)
            if !isnothing(ix)
                rr = direct_events[ix].rain_rate
                rain .+= (rr isa Matrix ? rr : fill(rr, grid_size))
            end
        end

        seq_ix = findlast(se -> se.timestamp <= t, seq)

        # Use dynamic equilibrium for directly leaking traps
        for trap in 1:n_traps
            leakage_state.leaking[trap] || continue
            haskey(drain_loc, trap) || continue
            leakage_state.leakage_start_time[trap] > t && continue

            V_leak = leakage_state.leakage_volume[trap]
            V_res = V_leak * sat
            dr = drain_rates[trap]

            inflow = _get_trap_inflow(trap, seq_ix)

            # Advance volume from previous interval.
            # Skip if the previous timestamp was before leakage started.
            if ti > 1 && sorted_ts[ti-1] >= leakage_state.leakage_start_time[trap]
                dt = t - sorted_ts[ti-1]
                prev_q = get(prev_inflow, trap, 0.0)
                if dt > 0
                    if prev_q > 0 && trap_volume[trap] < V_leak
                        # Was refilling
                        trap_volume[trap] = min(trap_volume[trap] + prev_q * dt, V_leak)
                    elseif prev_q <= 0 && trap_volume[trap] > V_res && has_finite_drainage
                        # Was draining
                        trap_volume[trap] = max(trap_volume[trap] - dr * dt, V_res)
                    end
                end
            end

            prev_inflow[trap] = inflow

            # Determine leakage rate based on current state
            if inflow > 0 && trap_volume[trap] >= V_leak
                # Equilibrium: column maintained, all inflow passes through
                rate = inflow * unit_factor
            elseif inflow > 0 && trap_volume[trap] < V_leak
                # Refilling: inflow goes to rebuilding column, excess passes through
                # refill_deficit = V_leak - trap_volume[trap]
                rate = 0.0
            elseif inflow <= 0 && trap_volume[trap] > V_res && has_finite_drainage
                # No inflow, draining: buoyancy drives CO2 through caprock
                rate = dr * unit_factor
            else
                # Fully drained or no activity
                rate = 0.0
            end

            rate > 0.0 && spread_rate!(rain, Tuple(drain_loc[trap]), rate, leakage_radius; regions=target_regions)
        end

        # Non-leaking draining traps has a constant-rate drainage
        if has_finite_drainage
            for trap in 1:n_traps
                leakage_state.leaking[trap] && continue  # handled above
                drain_rates[trap] == 0.0 && continue
                t0 = leakage_state.leakage_start_time[trap]
                (t >= t0 && t < t0 + T_res) || continue
                spread_rate!(rain, Tuple(drain_loc[trap]), drain_rates[trap] * unit_factor, leakage_radius; regions=target_regions)
            end
        end

        # Emit only when the combined rate actually changes.
        if isnothing(last_rain) || rain != last_rain
            push!(result, WeatherEvent(t, copy(rain)))
            last_rain = copy(rain)
        end
    end

    # Just for safety, if the last emitted event still carries a positive rate,
    # append a zero-rate event at the latest drainage end time.
    if !isempty(result) && any(result[end].rain_rate .> 0.0)
        max_t_end = -Inf
        for trap in 1:n_traps
            leakage_state.draining[trap] || continue
            t0 = leakage_state.leakage_start_time[trap]
            isfinite(t0) || continue
            if leakage_state.leaking[trap]
                # For leaking traps, drainage timing depends on inflow history.
                # Use cumulative_no_inflow_time + any remaining no-inflow period.
                # Conservative estimate: assume drainage could take up to T_res
                # from the last known state change.
                t_last_change = leakage_state.time_of_last_state_change[trap]
                if isfinite(t_last_change)
                    max_t_end = max(max_t_end, t_last_change + T_res)
                else
                    max_t_end = max(max_t_end, t0 + T_res)
                end
            else
                max_t_end = max(max_t_end, t0 + T_res)
            end
        end
        if isfinite(max_t_end) && max_t_end > result[end].timestamp
            push!(result, WeatherEvent(max_t_end, zeros(Float64, grid_size)))
        end
    end

    return result
end
