import Interpolations
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, subtrapsof, FilledAmount
using Distributions: Normal, LogNormal, Uniform, truncated

export volume_to_height, compute_leakage_volume
export find_leakage_location
export initialize_leakage_state

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
        leakage_volumes[trap_id] = leakage_volume === nothing ? 0.0 : leakage_volume
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