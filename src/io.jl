using NPZ
using Statistics
using CairoMakie

struct SleipnerTopography
    surfaces::Dict{String, Any}
    top_caprock::Array{Float64,2}
    sand_layers::Vector{Dict{String, Any}}
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    depth_min::Float64
    depth_max::Float64
end

function load_sleipner_topography(path::String = "sleipner/depth_surfaces/")
    println("\nLoading Sleipner depth surfaces...")

    # Load individual .npy files instead of .npz
    function load_surface(name::String)
        return npzread(joinpath(path, "$(name).npy"))
    end

    top_caprock = load_surface("Top_Caprock")
    top_sw = load_surface("TopSW")
    top_utsira = load_surface("TopUtsiraFm")
    base_utsira = load_surface("BaseUtsiraFm")
    thick_shale = load_surface("ThickShale")

    # Store all surfaces in a dictionary for compatibility
    surfaces = Dict{String, Any}(
        "Top_Caprock" => top_caprock,
        "TopSW" => top_sw,
        "TopUtsiraFm" => top_utsira,
        "BaseUtsiraFm" => base_utsira,
        "ThickShale" => thick_shale
    )

    sand_layers = []

    # L9: Shallowest sand (above thick shale)
    push!(sand_layers, Dict(
        "id" => 9,
        "name" => "L9",
        "top" => top_sw,
        "base" => thick_shale,
    ))

    # Load reflector surfaces
    for i in 1:7
        surfaces["Reflector$(i)"] = load_surface("Reflector$(i)")
        surfaces["Base_Reflector$(i)"] = load_surface("Base_Reflector$(i)")
    end

    # L8: First sand below thick shale
    push!(sand_layers, Dict(
        "id" => 8,
        "name" => "L8",
        "top" => top_utsira,
        "base" => surfaces["Reflector7"],
    ))

    # L7-L2: Sand layers between thin shales
    for i in 7:-1:2
        layer_num = i
        push!(sand_layers, Dict(
            "id" => layer_num,
            "name" => "L$(layer_num)",
            "top" => surfaces["Base_Reflector$(i)"],
            "base" => surfaces["Reflector$(i-1)"],
        ))
    end

    # L1: Deepest sand (below all thin shales)
    push!(sand_layers, Dict(
        "id" => 1,
        "name" => "L1",
        "top" => surfaces["Base_Reflector1"],
        "base" => base_utsira,
    ))

    # Sort by depth (shallowest to deepest by top surface mean)
    sand_layers = sort(sand_layers, by = x -> mean(x["top"]))

    nx, ny = size(top_caprock)
    depth_min = minimum(top_caprock)
    depth_max = maximum(base_utsira)

    # Bit hacky, but matches original grid spacing
    dx = 3200.0 / nx 
    dy = 5900.0 / ny

    return SleipnerTopography(
        surfaces,
        top_caprock,
        sand_layers,
        nx,
        ny,
        dx,
        dy,
        depth_min,
        depth_max,
    )
end

function create_domain_from_topography(topography::SleipnerTopography, dz::Float64)::Domain3D
    nx = topography.nx
    ny = topography.ny
    nz = Int(ceil((topography.depth_max - topography.depth_min) / dz))

    length_x = nx * topography.dx
    length_y = ny * topography.dy
    length_z = nz * dz

    Domain3D(
        nx,
        ny,
        nz,
        length_x,
        length_y,
        length_z,
        topography.depth_min,
        topography.depth_max,
    )
end
