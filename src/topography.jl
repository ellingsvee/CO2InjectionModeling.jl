"""
    AbstractTopography

Abstract base type for reservoir topography descriptions.

Subtypes must implement:
- [`get_sand_layers`](@ref)
- [`get_grid_dimensions`](@ref)
- [`get_grid_spacing`](@ref)
- [`get_depth_range`](@ref)

Default implementations are provided for [`get_caprock_surface`](@ref),
[`get_coordinate_origin`](@ref), and [`get_num_layers`](@ref).

See [`GenericTopography`](@ref) for the standard concrete implementation.
"""
abstract type AbstractTopography end

export AbstractTopography, GenericTopography
export get_sand_layers, get_grid_dimensions, get_grid_spacing, get_depth_range
export get_caprock_surface, get_coordinate_origin, get_num_layers
export create_domain

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

"""
    GenericTopography <: AbstractTopography

A topography constructed from raw surface arrays, suitable for use from R
or any context where dataset-specific loaders are not available.

# Fields
- `sand_layers`: Vector of layer Dicts with "name", "top", "base" keys
- `nx`, `ny`: Grid dimensions
- `dx`, `dy`: Grid spacing in meters
- `depth_min`, `depth_max`: Depth range in meters
- `caprock_surface`: Optional caprock top surface (nx × ny)
- `coordinate_origin`: (x, y) coordinate origin tuple
"""
struct GenericTopography <: AbstractTopography
    sand_layers::Vector{Dict{String,Any}}
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    depth_min::Float64
    depth_max::Float64
    caprock_surface::Union{Matrix{Float64},Nothing}
    coordinate_origin::Tuple{Float64,Float64}
end

# Convenience constructor with keyword arguments for optional fields
function GenericTopography(
    sand_layers::Vector{Dict{String,Any}},
    nx::Int, ny::Int,
    dx::Float64, dy::Float64,
    depth_min::Float64, depth_max::Float64;
    caprock_surface::Union{Matrix{Float64},Nothing}=nothing,
    coordinate_origin::Tuple{Float64,Float64}=(0.0, 0.0)
)
    GenericTopography(sand_layers, nx, ny, dx, dy, depth_min, depth_max, caprock_surface, coordinate_origin)
end

get_sand_layers(t::GenericTopography) = t.sand_layers
get_grid_dimensions(t::GenericTopography) = (t.nx, t.ny)
get_grid_spacing(t::GenericTopography) = (t.dx, t.dy)
get_depth_range(t::GenericTopography) = (t.depth_min, t.depth_max)
get_caprock_surface(t::GenericTopography) = t.caprock_surface
get_coordinate_origin(t::GenericTopography) = t.coordinate_origin
