"""
    AbstractTopography

Abstract base type for reservoir topography definitions.

Users implementing custom reservoir topographies should subtype this and implement
the required interface methods:
- `get_sand_layers(t)` - Returns Vector of layer Dicts with "name", "top", "base" keys
- `get_grid_dimensions(t)` - Returns (nx, ny) tuple
- `get_grid_spacing(t)` - Returns (dx, dy) tuple
- `get_depth_range(t)` - Returns (depth_min, depth_max) tuple

Optional methods with defaults:
- `get_caprock_surface(t)` - Returns caprock top surface or nothing
- `get_coordinate_origin(t)` - Returns (x, y) coordinate origin tuple
- `get_num_layers(t)` - Returns number of sand layers

Example:
```julia
struct MyTopography <: AbstractTopography
    layers::Vector{Dict{String,Any}}
    nx::Int; ny::Int; dx::Float64; dy::Float64
    depth_min::Float64; depth_max::Float64
end

CO2BatchFill.get_sand_layers(t::MyTopography) = t.layers
CO2BatchFill.get_grid_dimensions(t::MyTopography) = (t.nx, t.ny)
CO2BatchFill.get_grid_spacing(t::MyTopography) = (t.dx, t.dy)
CO2BatchFill.get_depth_range(t::MyTopography) = (t.depth_min, t.depth_max)

# Then use generic functions
domain = create_domain(my_topo, 1.0)
layers = analyze_base_surfaces(my_topo; boundary_condition=:closed)
```
"""
abstract type AbstractTopography end

export AbstractTopography
export get_sand_layers, get_grid_dimensions, get_grid_spacing, get_depth_range
export get_caprock_surface, get_coordinate_origin, get_num_layers
export create_domain

# =============================================================================
# Required interface methods (must be implemented by subtypes)
# =============================================================================

"""
    get_sand_layers(topography::AbstractTopography) -> Vector{Dict{String,Any}}

Return the sand layers of the reservoir. Each layer Dict must have:
- "name"::String - Layer name (e.g., "L1", "L2")
- "top"::Matrix{Float64} - Top surface depths (nx × ny)
- "base"::Matrix{Float64} - Base surface depths (nx × ny)
"""
function get_sand_layers end

"""
    get_grid_dimensions(topography::AbstractTopography) -> Tuple{Int, Int}

Return the grid dimensions (nx, ny) of the topography.
"""
function get_grid_dimensions end

"""
    get_grid_spacing(topography::AbstractTopography) -> Tuple{Float64, Float64}

Return the grid cell sizes (dx, dy) in meters.
"""
function get_grid_spacing end

"""
    get_depth_range(topography::AbstractTopography) -> Tuple{Float64, Float64}

Return the depth range (depth_min, depth_max) in meters.
depth_min is the shallowest depth, depth_max is the deepest.
"""
function get_depth_range end

# =============================================================================
# Optional interface methods (with default implementations)
# =============================================================================

"""
    get_caprock_surface(topography::AbstractTopography) -> Union{Matrix{Float64}, Nothing}

Return the caprock top surface depths, or nothing if not defined.
Default implementation returns nothing.
"""
get_caprock_surface(::AbstractTopography) = nothing

"""
    get_coordinate_origin(topography::AbstractTopography) -> Tuple{Float64, Float64}

Return the (x, y) coordinate origin for the grid (e.g., UTM coordinates).
Default implementation returns (0.0, 0.0).
"""
get_coordinate_origin(::AbstractTopography) = (0.0, 0.0)

"""
    get_num_layers(topography::AbstractTopography) -> Int

Return the number of sand layers.
Default implementation returns the length of get_sand_layers(topography).
"""
get_num_layers(t::AbstractTopography) = length(get_sand_layers(t))

# =============================================================================
# Generic domain creation
# =============================================================================

"""
    create_domain(topography::AbstractTopography, dz::Float64) -> Domain3D

Create a 3D domain from a topography definition.

# Arguments
- `topography`: Any AbstractTopography implementation
- `dz`: Vertical grid spacing in meters

# Returns
- `Domain3D`: A 3D domain structure suitable for simulation
"""
function create_domain(topography::AbstractTopography, dz::Float64)::Domain3D
    nx, ny = get_grid_dimensions(topography)
    dx, dy = get_grid_spacing(topography)
    depth_min, depth_max = get_depth_range(topography)

    nz = Int(ceil((depth_max - depth_min) / dz))

    length_x = nx * dx
    length_y = ny * dy
    length_z = nz * dz

    Domain3D(
        nx,
        ny,
        nz,
        length_x,
        length_y,
        length_z,
        depth_min,
        depth_max,
    )
end
