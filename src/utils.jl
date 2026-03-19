import Graphs
import Interpolations

export get_all_parents, get_all_descendants
export get_min_topography_elevation, get_trap_bottom_elevation
export volume_to_height, height_map
export get_grid_size, get_padded_grid_size, create_injection_rate, create_no_injection, volume_scale
export spread_rate!

"""
    get_all_parents(tstruct::TrapStructure, trap_id::Int)::Vector{Int}

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
    get_all_descendants(tstruct::TrapStructure, trap_id::Int)::Vector{Int}

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
    get_min_topography_elevation(trap_id::Int, tstruct::TrapStructure)::Float64

Get the minimum topography elevation in a trap's footprint.
This is the actual lowest point in the footprint, including the area covered by child traps.
"""
function get_min_topography_elevation(trap_id::Int, tstruct::TrapStructure)::Float64
    footprint = tstruct.footprints[trap_id]
    return minimum(tstruct.topography[footprint])
end

"""
    get_trap_bottom_elevation(trap_id::Int, tstruct::TrapStructure)::Float64

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
    volume_to_height(volume::Float64, trap_id::Int, z_vol_table::Tuple{Vector{Float64},Vector{Float64}}, tstruct::TrapStructure)::Float64

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

"""
    height_map(tstruct::TrapStructure, z_vol_tables::Vector{Tuple{Vector{Float64},Vector{Float64}}}, volumes::Vector{Float64};

Compute a 2-D CO2 column height map from per-trap volumes.

Optional keyword arguments support physically-correct height reporting for traps
that have leaked and are holding residually-trapped CO2.
"""
function height_map(tstruct, z_vol_tables, volumes;
    leaking::Union{Nothing,AbstractVector{Bool}}=nothing,
    leakage_heights::Union{Nothing,AbstractVector{Float64}}=nothing)
    num_traps = numtraps(tstruct)
    height_map = zeros(Float64, size(tstruct.regions))
    use_leakage_correction = !isnothing(leaking) && !isnothing(leakage_heights)
    for trap_id in 1:num_traps
        vol = volumes[trap_id]
        vol <= 0.0 && continue
        h = if use_leakage_correction && leaking[trap_id] && isfinite(leakage_heights[trap_id])
            leakage_heights[trap_id]
        else
            volume_to_height(vol, trap_id, z_vol_tables[trap_id], tstruct)
        end
        h <= 0.0 && continue
        min_topo = get_min_topography_elevation(trap_id, tstruct)
        water_level = min_topo + h
        for idx in tstruct.footprints[trap_id]
            cell_h = max(0.0, water_level - tstruct.topography[idx])
            height_map[idx] = max(height_map[idx], cell_h)
        end
    end
    return height_map
end

"""
    get_grid_size(layers) -> Tuple{Int,Int}

Return the (unpadded) grid size `(nx, ny)` from the first layer's trap structure topography.
Eliminates the need for manual padding arithmetic when setting up injection rates.
"""
function get_grid_size(layers::Vector{Layer})::Tuple{Int,Int}
    pad = layers[1].pad_width
    return size(layers[1].trap_structure.topography[1+pad:end-pad, 1+pad:end-pad])
end

"""
    get_grid_size(layer) -> Tuple{Int,Int}

Return the (unpadded) grid size `(nx, ny)` from a layer's trap structure topography.
Eliminates the need for manual padding arithmetic when setting up injection rates.
"""
function get_grid_size(layer::Layer)::Tuple{Int,Int}
    pad = layer.pad_width
    return size(layer.trap_structure.topography[1+pad:end-pad, 1+pad:end-pad])
end

"""
    get_padded_grid_size(layers) -> Tuple{Int,Int}

Return the (padded) grid size `(nx, ny)` from the first layer's trap structure topography.
Eliminates the need for manual padding arithmetic when setting up injection rates.
"""
function get_padded_grid_size(layers::Vector{Layer})::Tuple{Int,Int}
    return size(layers[1].trap_structure.topography)
end

"""
    get_padded_grid_size(layer) -> Tuple{Int,Int}

Return the (padded) grid size `(nx, ny)` from a layer's trap structure topography.
"""
function get_padded_grid_size(layer::Layer)::Tuple{Int,Int}
    return size(layer.trap_structure.topography)
end

"""
    spread_rate!(matrix, location, rate, radius)

Add `rate` to `matrix` spread evenly across all cells within Euclidean distance
`radius` of `location`. If `radius == 0`, the rate is added to the single cell.
"""
function spread_rate!(matrix::Matrix{Float64}, location::Tuple{Int,Int}, rate::Float64, radius::Int)
    if radius == 0
        matrix[location[1], location[2]] += rate
        return matrix
    end
    rows, cols = size(matrix)
    ci, cj = location
    # Collect all cells within the disc
    cells = Tuple{Int,Int}[]
    for i in max(1, ci - radius):min(rows, ci + radius)
        for j in max(1, cj - radius):min(cols, cj + radius)
            if (i - ci)^2 + (j - cj)^2 <= radius^2
                push!(cells, (i, j))
            end
        end
    end
    per_cell = rate / length(cells)
    for (i, j) in cells
        matrix[i, j] += per_cell
    end
    return matrix
end

"""
    create_injection_rate(layers, location_idx, rate) -> Matrix{Float64}

Create an injection rate matrix of the correct (padded) grid size with `rate` placed
at the specified *unpadded* `(i, j)` index. The padding offset for closed boundary
conditions is applied automatically.

For a zero-rate matrix, pass any index with `rate = 0.0`.

# Arguments
- `layers`: Vector of `Layer` from `analyze_base_surfaces`
- `location_idx`: `(i, j)` tuple in the original (unpadded) grid coordinates
- `rate`: Injection rate value to place at the specified location
"""
function create_injection_rate(layers::Vector{Layer}, location_idx::Tuple{Int,Int}, rate::Real; radius::Int=0)::Matrix{Float64}
    grid_size = get_padded_grid_size(layers)
    rate_matrix = zeros(grid_size)
    if rate != 0.0
        pad = layers[1].pad_width # 0 if :open, > 0 if :closed
        padded_idx = (location_idx[1] + pad, location_idx[2] + pad)
        spread_rate!(rate_matrix, padded_idx, Float64(rate), radius)
    end
    return rate_matrix
end

"""
    create_no_injection(layers) -> Vector{InjectionEvent}

Return `[InjectionEvent(0.0, zeros(grid_size))]` — the standard "no direct injection"
schedule for a non-injecting layer.
"""
function create_no_injection(layers::Vector{Layer})::Vector{InjectionEvent}
    grid_size = get_padded_grid_size(layers)
    return [InjectionEvent(0.0, zeros(grid_size))]
end

"""
    volume_scale(rp, domain) -> Float64

Return the conversion factor from SWIM volume units to physical m³.
Equivalent to `swim_volume_to_physical_volume(1.0, rp, domain)`.
"""
function volume_scale(rp::ReservoirProperties, domain::Domain3D)::Float64
    return swim_volume_to_physical_volume(1.0, rp, domain)
end