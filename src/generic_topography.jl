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

# Example
```julia
layers = [
    Dict{String,Any}("name" => "L1", "top" => top1, "base" => base1),
    Dict{String,Any}("name" => "L2", "top" => top2, "base" => base2),
]
topo = GenericTopography(layers, nx, ny, dx, dy, depth_min, depth_max)
```
"""
struct GenericTopography <: AbstractTopography
    sand_layers::Vector{Dict{String,Any}}
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    depth_min::Float64
    depth_max::Float64
    caprock_surface::Union{Matrix{Float64}, Nothing}
    coordinate_origin::Tuple{Float64, Float64}
end

# Convenience constructor without optional fields
function GenericTopography(
    sand_layers::Vector{Dict{String,Any}},
    nx::Int, ny::Int,
    dx::Float64, dy::Float64,
    depth_min::Float64, depth_max::Float64;
    caprock_surface::Union{Matrix{Float64}, Nothing}=nothing,
    coordinate_origin::Tuple{Float64, Float64}=(0.0, 0.0)
)
    GenericTopography(sand_layers, nx, ny, dx, dy, depth_min, depth_max, caprock_surface, coordinate_origin)
end

export GenericTopography

# AbstractTopography interface implementation
get_sand_layers(t::GenericTopography) = t.sand_layers
get_grid_dimensions(t::GenericTopography) = (t.nx, t.ny)
get_grid_spacing(t::GenericTopography) = (t.dx, t.dy)
get_depth_range(t::GenericTopography) = (t.depth_min, t.depth_max)
get_caprock_surface(t::GenericTopography) = t.caprock_surface
get_coordinate_origin(t::GenericTopography) = t.coordinate_origin
