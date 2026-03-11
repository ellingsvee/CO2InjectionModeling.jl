```@meta
EditURL = "../../../examples/sleipner_analysis.jl"
```

# Sleipner Analysis
In this example, we apply `CO2BatchFill` to analyse Sleipner, the first commercial-scale $\text{CO}_2$ storage project. We use publicly available depth surfaces downloaded from the Sleipner 2019 Benchmark Model.
Be aware that this example is intended for demonstration purposes only, and is not intended to match the true historical evolution of the Sleipner plume.

## Loading necessary packages

````@example sleipner_analysis
using CO2BatchFill
using SurfaceWaterIntegratedModeling
using CairoMakie # for visualization
using LaTeXStrings
using Random
using Statistics
using NPZ # for loading depth surfaces and feeder locations from NumPy files
````

## Setting up the layers and domain from the depth surfaces
### Loading data

````@example sleipner_analysis
data_dir = CO2BatchFill.datapath_testdata()
depth_surfaces = joinpath(data_dir, "sleipner", "depth_surfaces")
````

### Defining the topography struct
The `SleipnerTopography` struct implements the `AbstractTopography` interface required by `CO2BatchFill`.

````@example sleipner_analysis
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
CO2BatchFill.get_sand_layers(t::SleipnerTopography) = t.sand_layers
CO2BatchFill.get_grid_dimensions(t::SleipnerTopography) = (t.nx, t.ny)
CO2BatchFill.get_grid_spacing(t::SleipnerTopography) = (t.dx, t.dy)
CO2BatchFill.get_depth_range(t::SleipnerTopography) = (t.depth_min, t.depth_max)
CO2BatchFill.get_caprock_surface(t::SleipnerTopography) = t.top_caprock
CO2BatchFill.get_coordinate_origin(::SleipnerTopography) = (436885.1562, 6469129.5) # UTM coordinates of the bottom-left corner of the grid
````

### Loading and processing the depth surfaces
This function loads the depth surfaces from the provided `.npy` files, constructs the sand layers based on the reflector and shale surfaces.

````@example sleipner_analysis
function load_sleipner_topography(path::String)
    # Load individual .npy files instead of .npz
    # Note: We reverse the Y-axis (2nd dimension) to match UTM coordinate convention
    # where northing increases upward, but array indices increase downward
    function load_surface(name::String)
        data = npzread(joinpath(path, "$(name).npy"))
        return reverse(data, dims=2)  # Flip Y-axis to correct orientation
    end

    top_caprock = load_surface("Top_Caprock")
    top_sw = load_surface("TopSW")
    top_utsira = load_surface("TopUtsiraFm")
    base_utsira = load_surface("BaseUtsiraFm")
    thick_shale = load_surface("ThickShale")

    # Store all surfaces in a dictionary for compatibility
    surfaces = Dict{String,Any}(
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
    sand_layers = sort(sand_layers, by=x -> mean(x["top"]))

    nx, ny = size(top_caprock)
    depth_min = minimum(top_caprock)
    depth_max = maximum(base_utsira)

    # Grid spacing from metadata cell-center extent: xmin=436910.1562 to xmax=440110.1562 (3200m),
    # ymin=6469154.5 to ymax=6475054.5 (5900m). These are center-to-center spans, so divide by (n-1).
    dx = 3200.0 / (nx - 1)
    dy = 5900.0 / (ny - 1)

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
````

Loading and analyzing the topography

````@example sleipner_analysis
topography = load_sleipner_topography(depth_surfaces)
layers = analyze_base_surfaces(topography; boundary_condition=:closed)
````

Create the `Domain3D` for the simulation.

````@example sleipner_analysis
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
domain = create_domain_from_topography(topography, 1.0)
````

## Injection events
Based on the historical injection data for Sleipner, we create a time series of injection events.

````@example sleipner_analysis
function generate_sleipner_injection_events(
    layers,
    injection_cell::CartesianIndex
)::Vector{Vector{InjectionEvent}}

    # Historical annual injection rates from 1996-2010 (Mt/year)
    # Source: Sleipner 2019 Benchmark model
    annual_rates_mt = [
        0.07,  # 1996
        0.67,  # 1997
        0.85,  # 1998
        0.94,  # 1999
        0.94,  # 2000
        1.02,  # 2001
        0.96,  # 2002
        0.92,  # 2003
        0.76,  # 2004
        0.87,  # 2005
        0.83,  # 2006
        0.93,  # 2007
        0.82,  # 2008
        0.86,  # 2009
        0.76   # 2010
    ]

    # Convert Mt/year to m^3/year
    co2_density_l1 = 460.0
    annual_rates_m3_per_year = annual_rates_mt .* 1e9 ./ co2_density_l1

    # Create injection events for bottom layer (L1)
    # Time points are cumulative: start at year 0, events mark end of each year
    n_events = length(annual_rates_mt)
    bottom_layer_events = Vector{InjectionEvent}(undef, n_events)

    # Get the grid size from the bottom layer
    grid_size = size(layers[1].trap_structure.topography)

    for (i, rate) in enumerate(annual_rates_m3_per_year)
        # Time in years (0, 1, 2, ..., 14)
        time = float(i - 1)

        # Create injection rate field (only inject at specified cell)
        injection_rate = zeros(grid_size)
        injection_rate[injection_cell] = rate
        bottom_layer_events[i] = InjectionEvent(time, injection_rate)
    end

    # Create zero injection events for all other layers
    zero_injection = zeros(grid_size)
    zero_event = [InjectionEvent(0.0, zero_injection)]

    # Assemble injection events for all layers
    n_layers = length(layers)
    injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)
    injection_events[1] = bottom_layer_events
    for i in 2:n_layers
        injection_events[i] = zero_event
    end

    return injection_events
end

injection_events = generate_sleipner_injection_events(
    layers,
    CartesianIndex(div(topography.nx, 2), div(topography.ny, 2)) # Inject at center of domain
)
````

## Reservoir properties
We choose reservoir properties based on typical values for the Utsira Formation, which is the main storage formation at Sleipner.
For simplicity, we assume the same properties for all layers, except for an impermeable top layer.

````@example sleipner_analysis
function generate_reservoir_properties_for_sleipner_layers()::Vector{ReservoirProperties}
    n_layers = 9

    # Common reservoir properties for all layers
    sand_porosity::Float64 = 0.4
    sand_residual_co2_saturation::Float64 = 0.2
    sand_irreducible_water_saturation::Float64 = 0.3
    shale_pressure_threshold::Float64 = 98000.0
    residual_leakage_time::Float64 = 5.0 # years

    shale_pressure_thresholds = fill(shale_pressure_threshold, n_layers)
    shale_pressure_thresholds[end] = Inf  # Top layer impermeable

    # Density values from L1 up to L9 (from paper)
    brine_density = 1020.0
    co2_density = 460.0  # Average CO2 density

    # Create ReservoirProperties for each layer
    reservoir_properties = Vector{ReservoirProperties}(undef, n_layers)
    for i in 1:n_layers
        reservoir_properties[i] = ReservoirProperties(
            sand_porosity,
            sand_residual_co2_saturation,
            sand_irreducible_water_saturation,
            shale_pressure_thresholds[i],
            residual_leakage_time;
            # leakage_height computed automatically from pressure
            brine_density=brine_density,
            co2_density=co2_density
        )
    end
    return reservoir_properties
end

all_layers_rp = generate_reservoir_properties_for_sleipner_layers()
````

## Run simulation
Finally, we can run the simulation using the `fill_layers` function

````@example sleipner_analysis
seqs, leakage_states, weather_events_per_layer = fill_layers(
    layers, domain, all_layers_rp, injection_events; verbose=false)
````

Based on the simulation results, we generate yearly snapshots

````@example sleipner_analysis
t_end = 15.0
timepoints = collect(range(0.0, stop=t_end, length=30))
multi_snaps = generate_multi_layer_snapshots(
    layers, seqs, leakage_states, weather_events_per_layer, timepoints)
````

Print summary at final snapshot. See that the $\text{CO}_2$ has migrated up to the top layer.

````@example sleipner_analysis
print_summary(multi_snaps[end])
````

## Visualization
Optional defaults for making the plots look nicer.

````@example sleipner_analysis
update_theme!(
    merge(
        theme_latexfonts(),
        Theme(fontsize=25)
    )
)
````

Plot the plume extents at the time of the final snapshot.

````@example sleipner_analysis
injection_location_loc = (div(topography.nx, 2) * domain.dx, div(topography.ny, 2) * domain.dy)
plot_multi_layer(
    layers, multi_snaps[end], domain;
    max_co2_height=ceil(round(all_layers_rp[1].leakage_height, digits=2)),
    show_contours=true,
    show_labels=true,
    contour_levels=20,
    major_contour_every=5,
    contour_opacity=1.0,
    figure_size=(500 * 3, 600 * 3),
    figure_layout=(3, 3),
    colormap=:Blues,
    injection_locations=[injection_location_loc],
    show_leakage_locations=true,
    show_extents=true,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

