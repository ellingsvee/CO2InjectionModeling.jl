using NPZ
using Statistics

export SleipnerTopography

"""
    SleipnerTopography <: AbstractTopography

Topography definition for the Sleipner CO2 storage site in the North Sea.

Contains depth surfaces for the Utsira Formation sand layers and caprock,
loaded from NPZ files.

# Fields
- `surfaces`: Dict containing all loaded surface arrays
- `top_caprock`: Top of caprock surface (nx × ny)
- `sand_layers`: Vector of layer Dicts with "name", "top", "base" keys
- `nx`, `ny`: Grid dimensions
- `dx`, `dy`: Grid spacing in meters
- `depth_min`, `depth_max`: Depth range in meters
"""
struct SleipnerTopography <: AbstractTopography
    surfaces::Dict{String,Any}
    top_caprock::Array{Float64,2}
    sand_layers::Vector{Dict{String,Any}}
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    depth_min::Float64
    depth_max::Float64
end

# =============================================================================
# AbstractTopography interface implementation
# =============================================================================

CO2BatchFill.get_sand_layers(t::SleipnerTopography) = t.sand_layers
CO2BatchFill.get_grid_dimensions(t::SleipnerTopography) = (t.nx, t.ny)
CO2BatchFill.get_grid_spacing(t::SleipnerTopography) = (t.dx, t.dy)
CO2BatchFill.get_depth_range(t::SleipnerTopography) = (t.depth_min, t.depth_max)
CO2BatchFill.get_caprock_surface(t::SleipnerTopography) = t.top_caprock

# Sleipner-specific coordinate origin (UTM Zone 31N)
CO2BatchFill.get_coordinate_origin(::SleipnerTopography) = (436800.0, 6468100.0)
