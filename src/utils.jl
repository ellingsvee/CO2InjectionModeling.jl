import Graphs
import Interpolations

export get_all_parents, get_all_descendants
export get_min_topography_elevation, get_trap_bottom_elevation
export volume_to_height, height_map

"""
Get all parent trap IDs for a given trap ID in a trap structure.
"""
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

"""
Get the minimum topography elevation in a trap's footprint.
This is the actual lowest point in the footprint, including the area covered by child traps.
"""
function get_min_topography_elevation(trap_id::Int, tstruct::TrapStructure)::Float64
    footprint = tstruct.footprints[trap_id]
    return minimum(tstruct.topography[footprint])
end

"""
Get the effective bottom elevation of a trap for volume interpolation purposes.

NOTE: This is used for z_vol_table interpolation, NOT for CO2 height calculation.
For height calculation (leakage detection), use get_min_topography_elevation instead.
"""
function get_trap_bottom_elevation(trap_id::Int, tstruct::TrapStructure)::Float64
    trap_bottom_elevation = get_min_topography_elevation(trap_id, tstruct)

    # For parent traps, the effective bottom is above child spillpoints
    children = subtrapsof(tstruct, trap_id)
    if !isempty(children)
        child_spillpoint_elev = tstruct.spillpoints[children[1]].elevation
        trap_bottom_elevation = max(trap_bottom_elevation, child_spillpoint_elev)
    end

    return trap_bottom_elevation
end




"""
Convert a volume in a trap to the CO2 column height above the topography bottom.

For parent traps with filled children, the CO2 column extends from the true
topography minimum (in the child footprints) through the children up to the
current water level.
"""
function volume_to_height(
    volume::Float64,
    trap_id::Int,
    z_vol_table::Tuple{Vector{Float64},Vector{Float64}},
    tstruct::TrapStructure
)::Float64

    min_topography_elevation = get_min_topography_elevation(trap_id, tstruct)

    zvals, vvals = z_vol_table

    if length(zvals) == 1 || volume <= 0.0
        return 0.0
    end

    # Create interpolation function from volume to z
    v2z = Interpolations.linear_interpolation(vvals, zvals, extrapolation_bc=Interpolations.Line())

    water_level = v2z(volume)
    height = max(0.0, water_level - min_topography_elevation)

    return height
end

function height_map(tstruct, z_vol_tables, volumes)
    num_traps = numtraps(tstruct)
    height_map = zeros(Float64, size(tstruct.regions))
    for trap_id in 1:num_traps
        vol = volumes[trap_id]
        vol <= 0.0 && continue
        h = volume_to_height(vol, trap_id, z_vol_tables[trap_id], tstruct)
        h <= 0.0 && continue
        for idx in tstruct.footprints[trap_id]
            height_map[idx] = max(height_map[idx], h)
        end
    end
    return height_map
end