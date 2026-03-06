using SurfaceWaterIntegratedModeling
export Layer, add_boundary_wall, analyze_base_surfaces

"""
    Layer

Represents one geological storage layer after trap analysis.

# Fields
- `name`: Layer name (e.g. `"Storage layer 1"`)
- `trap_structure`: SWIM `TrapStructure` produced by `spillanalysis`
- `boundary_condition`: `:open` (CO2 can spill out of the domain) or
  `:closed` (boundary walls added via padding)

Constructed automatically by [`analyze_base_surfaces`](@ref).
"""
struct Layer
    name::String
    trap_structure::TrapStructure
    boundary_condition::Symbol  # :open or :closed
end

"""
    add_boundary_wall(surface, pad_width) -> Matrix{Float64}

Surround a topography surface with a wall of very high elevation so that CO2
cannot spill out of the domain (closed boundary condition).

A ring of `pad_width` cells is added on every side, set to
`maximum(surface) + 1e6` m, ensuring that no spill path crosses the boundary.

# Arguments
- `surface`: `(nx, ny)` depth matrix (positive downward)
- `pad_width`: Number of cells to add on each side

# Returns
- Padded `(nx + 2*pad_width, ny + 2*pad_width)` matrix
"""
function add_boundary_wall(surface::Matrix{<:Real}, pad_width::Int)
    nx, ny = size(surface)
    padded_nx = nx + 2 * pad_width
    padded_ny = ny + 2 * pad_width
    padded_surface = zeros(Float64, padded_nx, padded_ny)
    padded_surface[pad_width+1:pad_width+nx, pad_width+1:pad_width+ny] .= surface

    # Create walls by setting boundary cells to very high elevation
    wall_elevation = maximum(surface) + 1e6

    # Top, bottom, left and right boundaries
    padded_surface[1:pad_width, :] .= wall_elevation
    padded_surface[end-pad_width+1:end, :] .= wall_elevation
    padded_surface[:, 1:pad_width] .= wall_elevation
    padded_surface[:, end-pad_width+1:end] .= wall_elevation

    return padded_surface
end

"""
    analyze_base_surfaces(topography; boundary_condition=:open, pad_width=2) -> Vector{Layer}

Analyze all sand layers in `topography` using SWIM's spill-analysis algorithm
and return a `Vector{Layer}` ordered **deepest first** (index 1 = injection layer).

Layers are provided shallowest-first in [`GenericTopography`](@ref) but are
reversed here so that `layers[1]` corresponds to the deepest (injection) layer
and `fill_layers` can propagate leakage upward in index order.

# Arguments
- `topography`: Any [`AbstractTopography`](@ref) (e.g. [`GenericTopography`](@ref))
- `boundary_condition`: `:open` — CO2 can leave the domain; `:closed` — boundary
  walls are added via [`add_boundary_wall`](@ref) with `pad_width` cells
- `pad_width`: Width of the boundary wall in grid cells (only used when
  `boundary_condition = :closed`)

# Returns
- `Vector{Layer}`, index 1 = deepest / injection layer

# Example
```julia
layers = analyze_base_surfaces(topo; boundary_condition=:closed, pad_width=2)
```
"""
function analyze_base_surfaces(topography::AbstractTopography; boundary_condition::Symbol=:open, pad_width::Int=2)::Vector{Layer}
    @assert boundary_condition in [:open, :closed] "boundary_condition must be :open or :closed"


    # Get topography properties via interface methods
    sand_layers = get_sand_layers(topography)
    nx, ny = get_grid_dimensions(topography)
    dx, dy = get_grid_spacing(topography)


    # The spillanalysis function in SWIM expects the physical dimensions of the domain.
    pad = boundary_condition == :closed ? pad_width : 0
    length_x = nx * dx + 2 * pad * dx
    length_y = ny * dy + 2 * pad * dy

    # Initialize empty vector (will grow as we push)
    layers = Vector{Layer}()

    # Iterate over each sand layer to create Layer structs
    for layer in reverse(sand_layers)
        layer_name = layer["name"]
        surface = layer["top"]

        # Add boundary walls if closed BC
        if boundary_condition == :closed
            surface = add_boundary_wall(surface, pad_width)
        end

        # Use spillanalysis from SWIM
        trap_structure = spillanalysis(surface, lengths=(length_x, length_y))
        push!(layers, Layer(layer_name, trap_structure, boundary_condition))
    end

    return layers
end