using SurfaceWaterIntegratedModeling
export Layer, add_boundary_wall, analyze_base_surfaces

struct Layer
    name::String
    trap_structure::TrapStructure
    boundary_condition::Symbol  # :open or :closed
end

"""
Add wall padding around a topography surface for closed boundary conditions.
"""
function add_boundary_wall(surface::Matrix{<:Real})
    pad_width = 1  # For closed BC, we add a 1-cell wall around the edges    

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
Use SWIM to analyze the base surfaces
"""
function analyze_base_surfaces(topography::AbstractTopography; boundary_condition::Symbol=:open)::Vector{Layer}
    @assert boundary_condition in [:open, :closed] "boundary_condition must be :open or :closed"

    # Determine padding width: 0 for open BC, 1 for closed BC

    # Get topography properties via interface methods
    sand_layers = get_sand_layers(topography)
    nx, ny = get_grid_dimensions(topography)
    dx, dy = get_grid_spacing(topography)


    # The spillanalysis function in SWIM expects the physical dimensions of the domain.
    length_x = nx * dx + (boundary_condition == :closed ? 2 * dx : 0)
    length_y = ny * dy + (boundary_condition == :closed ? 2 * dy : 0)

    # Initialize empty vector (will grow as we push)
    layers = Vector{Layer}()

    # Iterate over each sand layer to create Layer structs
    for layer in reverse(sand_layers)
        layer_name = layer["name"]
        surface = layer["top"]

        # Add boundary walls if closed BC
        if boundary_condition == :closed
            surface = add_boundary_wall(surface)
        end

        # Use spillanalysis from SWIM
        trap_structure = spillanalysis(surface, lengths=(length_x, length_y))
        push!(layers, Layer(layer_name, trap_structure, boundary_condition))
    end

    return layers
end