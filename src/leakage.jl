import Interpolations
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, subtrapsof, FilledAmount
using Distributions: Normal, LogNormal, Uniform, truncated

export volume_to_height, compute_leakage_volume
export find_leakage_location, find_next_leakage_event
export initialize_leakage_state, create_empty_leakage_state
export LeakageTimeEstimate, get_initial_leakage_time_estimates, update_leakage_time_estimates!

"""
Stores the estimated time when a trap will reach its leakage volume. Similar to SWIM's ChangeTimeEstimate.
"""
struct LeakageTimeEstimate
    trap_id::Int
    min_time::Float64  # Earliest possible leakage time
    max_time::Float64  # Latest possible leakage time
end

"""
    find_leakage_location(trap_id::Int, tstruct::TrapStructure) -> CartesianIndex{2}

Find the grid cell where leakage occurs for a trap.
This is the lowest point (minimum topography elevation) in the trap footprint.
"""
function find_leakage_location(trap_id::Int, tstruct::TrapStructure)::CartesianIndex{2}
    footprint = tstruct.footprints[trap_id]
    topo_vals = tstruct.topography[footprint]
    min_idx = argmin(topo_vals)
    linear_idx = footprint[min_idx]
    return CartesianIndices(size(tstruct.topography))[linear_idx]
end


"""
Compute the volume at which a trap reaches the leakage height threshold.

The CO2 column height is measured from the TRUE topography bottom (including child traps)
to the current CO2 level. This accounts for CO2 that fills child traps
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

    # The lowest point in the trap footprint (including child traps)
    min_topography_elevation = get_min_topography_elevation(trap_id, tstruct)

    # Leakage occurs when level reaches min_topography_elevation + leakage_height
    leakage_elevation = min_topography_elevation + leakage_height

    # Get spillpoint elevation (maximum fill level)
    spillpoint_elevation = tstruct.spillpoints[trap_id].elevation

    # If leakage elevation is above spillpoint, trap spills before it can leak
    if leakage_elevation >= spillpoint_elevation
        return nothing
    end

    # Use z_vol_table to find volume at leakage_elevation
    zvals, vvals = z_vol_table

    if leakage_elevation <= zvals[1]
        # Leakage elevation is below the trap's z_vol_table minimum.
        # For parent traps, zvals[1] is the child spillpoint.
        # This means leakage would occur while children are still filling,
        # so this trap's own volume at leakage is 0.
        return 0.0
    end

    # Create interpolation function from z to volume
    z2v = Interpolations.linear_interpolation(zvals, vvals, extrapolation_bc=Interpolations.Line())

    return z2v(leakage_elevation)
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

    # For now...can eventually support per-trap leakage heights
    trap_leakage_heights = fill(reservoir_properties.leakage_height, num_traps)

    # Precompute leakage volumes for all traps using their specific heights
    leakage_volumes = zeros(Float64, num_traps)
    for trap_id in 1:num_traps
        trap_height = trap_leakage_heights[trap_id]
        leakage_volume = compute_leakage_volume(trap_id, z_vol_tables[trap_id], tstruct, trap_height)
        leakage_volumes[trap_id] = leakage_volume === nothing ? Inf : leakage_volume
    end

    return LeakageState(
        fill(false, num_traps),           # leaking
        fill(false, num_traps),           # draining
        leakage_volumes,                   # leakage_volume
        fill(Inf, num_traps),             # leakage_start_time
        LeakageRecord[],                   # leakage_records
        trap_leakage_heights,              # leakage_height 
        fill(0.0, num_traps),             # initial_volume_at_leak (0 until leakage starts)
        residual_saturation,               # residual_saturation
        residual_leakage_time              # residual_leakage_time
    )
end

"""
An empty leakage state for layers with no traps (e.g., non-dome layers).
"""
function create_empty_leakage_state()::LeakageState
    return LeakageState(
        Bool[],
        Bool[],
        Float64[],
        Float64[],
        LeakageRecord[],
        Float64[],
        Float64[],
        0.0,
        0.0
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
        _compute_leakage_time_estimate(trap_id, cur_amounts, cur_time, rateinfo,
            leakage_state, filled_traps, tstruct)
        for trap_id in 1:numtraps(tstruct)
    ]
end

"""
Compute when a trap will reach its leakage threshold.
"""
function _compute_leakage_time_estimate(
    trap_id::Int,
    cur_amounts::Vector{FilledAmount},
    cur_time::Float64,
    rateinfo,  # SWIM's RateInfo
    leakage_state::LeakageState,
    filled_traps::Vector{Bool},
    tstruct::TrapStructure
)::LeakageTimeEstimate
    # If trap has no leakage volume (e.g., spills before reaching leakage height), return Inf
    no_leakage_estimate = LeakageTimeEstimate(trap_id, Inf, Inf)

    # Already leaking - no future leakage event
    if leakage_state.leaking[trap_id]
        return no_leakage_estimate
    end

    # Get target volume for leakage
    target_vol = leakage_state.leakage_volume[trap_id]

    # If target is Inf, trap cannot leak (spills before reaching leakage height)
    if target_vol == Inf
        return no_leakage_estimate
    end

    # Get current volume
    current_vol = cur_amounts[trap_id].amount

    # If trap is not accumulating (children not filled), cannot estimate leakage
    children = subtrapsof(tstruct, trap_id)
    if !all(filled_traps[children])
        return no_leakage_estimate
    end

    # Check if already at or above target
    if target_vol == 0.0 && current_vol > 0.0
        # Trap has CO2 and should leak immediately (already above threshold)
        @warn "Trap $trap_id has leakage volume of 0 but already contains CO2. Treat it as leaking immediately."
        return LeakageTimeEstimate(trap_id, cur_time, cur_time)
    elseif current_vol >= target_vol
        return LeakageTimeEstimate(trap_id, cur_time, cur_time)
    end

    # TODO: Think this logic is correct. Assume we dont need to account for any outflow here.

    # Compute time using inflow rate
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
Find the next leakage event (minimum leakage time).
Returns (time, trap_id). Returns (Inf, 0) is no future leakage events.
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
        leakage_time_est[trap_id] = _compute_leakage_time_estimate(
            trap_id, cur_amounts, cur_time, rateinfo,
            leakage_state, filled_traps, tstruct
        )
    end
end

"""
Compute the current stored volume in a draining trap, accounting for residual drainage.

The volume decreases linearly from initial_volume to residual_volume over the
residual_leakage_time period after drainage starts.

A trap is "draining" if it's either:
1. Directly leaking (reached leakage threshold, edge=0), or
2. A filled ancestor whose CO2 flows through a leaking trap

Returns:
- For draining traps: Current volume accounting for drainage
- For non-draining traps: `nothing` (caller should use actual volume from SWIM)

If residual_leakage_time is 0 or Inf, no drainage occurs (returns initial volume).
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