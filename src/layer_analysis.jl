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
Use SWIM to analyze the base surfaces
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