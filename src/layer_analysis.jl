using SurfaceWaterIntegratedModeling

struct Layer
    name::String
    trap_structure::TrapStructure
end


function analyze_base_surfaces(topography::SleipnerTopography)::Vector{Layer}

    # Initialize empty vector (will grow as we push)
    layers = Vector{Layer}()

    # Iterate over each sand layer to create Layer structs
    for layer in reverse(topography.sand_layers)
        layer_name = layer["name"]
        base_surface = layer["top"]
        trap_structure = spillanalysis(base_surface, lengths = (topography.nx * topography.dx, topography.ny * topography.dy));
        push!(layers, Layer(layer_name, trap_structure))
    end

    return layers
end
