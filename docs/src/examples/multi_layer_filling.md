```@meta
EditURL = "../../../examples/multi_layer_filling.jl"
```

# Multi-layer filling

This example demonstrates the core functionality of `CO2BatchFill`. We generate a synthetic subsurface domain with three layers, simulate $\text{CO}_2$ injection into the deepest layer, and analyze the resulting $\text{CO}_2$ distribution and migration over time.

## Loading necessary packages

````@example multi_layer_filling
using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields # for generating synthetic surfaces
using CairoMakie # for visualization
using LaTeXStrings
using Random
````

## Selecting physical properties

````@example multi_layer_filling
const N_LAYERS = 3
const NX, NY = 100, 100
const LENGTH_X = 1000.0
const LENGTH_Y = 1000.0
const DX, DY = LENGTH_X / NX, LENGTH_Y / NY
const LAYER_THICK = 10.0     # m (vertical thickness of each sand layer)
const RESIDUAL_TRAPPING = 0.4 # fraction of CO2 volume that becomes immobile when a layer fills
const CAPILLARY_ENTRY_PRESSURE = 25_000.0
const boundary_condition = :closed
````

The `ReservoirProperties` struct defines the physical properties of each layer.

````@example multi_layer_filling
rp = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, CAPILLARY_ENTRY_PRESSURE, 5.0)
rp_caprock = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, Inf, 5.0) # Infinite entry pressure for caprock layer
all_layers_rp = [rp, rp, rp_caprock] # properties for each layer, in order from deepest to shallowest
````

Based on the `CAPILLARY_ENTRY_PRESSURE`, the `leakage_height` is the maximum $\text{CO}_2$ column height that can be supported before migration occurs.

````@example multi_layer_filling
println("Leakage height threshold: $(round(rp.leakage_height, digits=2)) m")
````

## Generating synthetic topography
We create three synthetic layers using Gaussian random fields (GRFs) with a Matern covariance function. This will represent the topographies of the shale-layers.

````@example multi_layer_filling
Random.seed!(101) # for reproducibility
cov = CovarianceFunction(2, Matern(200, 2, σ=3.0))
pts = range(0.0, stop=(NX - 1) * DX, step=DX)

function sample_surface(base_depth)
    sample(GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)) .+ base_depth
end

surf_L1 = sample_surface(950.0)   # deepest  (injection layer)
surf_L2 = sample_surface(900.0)
surf_L3 = sample_surface(850.0)
````

We create a `GenericTopography` from the three surfaces, which implements the `AbstractTopography` interface required by `CO2BatchFill`.

````@example multi_layer_filling
sand_layers = [
    Dict{String,Any}("name" => "Storage layer 3", "top" => surf_L3, "base" => surf_L3 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 2", "top" => surf_L2, "base" => surf_L2 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 1", "top" => surf_L1, "base" => surf_L1 .+ LAYER_THICK),
]
topo = GenericTopography(sand_layers, NX, NY, DX, DY, minimum(surf_L3), maximum(surf_L1) + LAYER_THICK)
domain = create_domain(topo, 0.1)
````

Using the `spillanalysis` algorithm from SWIM, the `analyze_base_surfaces` function identifies the spill points and leakage paths for each layer. The `boundary_condition` argument controls whether $\text{CO}_2$ can leave the domain (open BC) or if boundary walls are added to prevent outflow (closed BC).

````@example multi_layer_filling
layers = analyze_base_surfaces(topo; boundary_condition=boundary_condition)
````

## Defining injection schedule
We define an injection schedule where $\text{CO}_2$ is injected into the deepest layer (Layer 1) at a constant rate over ten years. The upper layers (Layers 2–3) have no direct injection. The injection schedule is defined as a matrix specifying the injection rate for each layer. The `create_injection_rate` and `create_no_injection` are utility functions that generate the appropriate injection rate matrices based on grid dimensions.

````@example multi_layer_filling
TOTAL_RATE = 80_000.0
INJECTION_END = 10.0 # injection up to 10 years, then stop
rate_L1 = create_injection_rate(layers, (div(NX, 2), div(NY, 2)), TOTAL_RATE)
null_rate = create_no_injection(layers)
````

We create an `injection_events` vector with one entry per layer, where each entry is a vector of `InjectionEvent`s.

````@example multi_layer_filling
injection_events = [
    [InjectionEvent(0.0, rate_L1), InjectionEvent(INJECTION_END, create_injection_rate(layers, (1, 1), 0.0))], # Layer 1
    null_rate, # Layer 2
    null_rate, # Layer 3
]
````

## Run simulation
The `fill_layers` function runs the simulation.

````@example multi_layer_filling
seqs, leakage_states, weather_events_per_layer = fill_layers(
    layers, domain, all_layers_rp, injection_events; verbose=false)
````

After the simulation is complete, we can generate snapshots at arbitrary time points. `generate_multi_layer_snapshots` creates a vector of `MultiLayerSnapshot`s containing information about the total $\text{CO}_2$ distribution across all layers, and a `LayerSnapshot`s for each individual layer.

````@example multi_layer_filling
t_end = 15.0
timepoints = collect(range(0.0, stop=t_end, length=30))
multi_snaps = generate_multi_layer_snapshots(
    layers, seqs, leakage_states, weather_events_per_layer, timepoints)
````

Print summary at final snapshot

````@example multi_layer_filling
print_summary(multi_snaps[end])
````

## Visualization
Optional defaults for making the plots look nicer.

````@example multi_layer_filling
update_theme!(
    merge(
        theme_latexfonts(),
        Theme(fontsize=25)
    )
)
````

First plot the $\text{CO}_2$ plumes at the time of the final snapshot. The `show_extents` argument controls whether to visualize the $\text{CO}_2$ heights, or simply the extent of the plume

````@example multi_layer_filling
injection_location_loc = (div(NX, 2) * DX, div(NY, 2) * DY)
plot_multi_layer(
    layers, multi_snaps[end], domain;
    max_co2_height=ceil(round(rp.leakage_height, digits=2)),
    show_contours=true,
    show_labels=true,
    contour_levels=20,
    major_contour_every=5,
    contour_opacity=1.0,
    figure_size=(500 * N_LAYERS, 500),
    colormap=:Blues,
    injection_locations=[injection_location_loc],
    show_leakage_locations=true,
    show_extents=false,
    colorbar_label="Column height",
)
````

Lastly, the `plot_multi_layer_volumes_timeseries` function plots the total $\text{CO}_2$ volume in each layer over time. The `vol_scale` argument can be used to convert the volumes to physical units (e.g. cubic meters) and adjust the y-axis scale accordingly.

````@example multi_layer_filling
_phys_scale = volume_scale(rp, domain)
plot_multi_layer_volumes_timeseries(
    multi_snaps;
    vol_scale=_phys_scale / 1e5,
    ylabel="Volume 10^5",
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

