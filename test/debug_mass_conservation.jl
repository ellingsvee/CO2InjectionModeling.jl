using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields # for generating synthetic surfaces
using LaTeXStrings
using Random


# ## Selecting physical properties
N_LAYERS = 3
NX, NY = 100, 100
LENGTH_X = 1000.0
LENGTH_Y = 1000.0
DX, DY = LENGTH_X / NX, LENGTH_Y / NY
LAYER_THICK = 10.0     # m (vertical thickness of each sand layer)
RESIDUAL_TRAPPING = 0.4 # fraction of CO2 volume that becomes immobile when a layer fills
CAPILLARY_ENTRY_PRESSURE = 25_000.0
boundary_condition = :closed

# The `ReservoirProperties` struct defines the physical properties of each layer.
rp = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, CAPILLARY_ENTRY_PRESSURE, 5.0)
rp_caprock = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, Inf, 5.0) # Infinite entry pressure for caprock layer 
all_layers_rp = [rp, rp, rp_caprock] # properties for each layer, in order from deepest to shallowest

# ## Generating synthetic topography
Random.seed!(101) # for reproducibility
cov = CovarianceFunction(2, Matern(50, 2, σ=3.0))
pts = range(0.0, stop=(NX - 1) * DX, step=DX)

function sample_surface(base_depth)
    sample(GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)) .+ base_depth
end

surf_L1 = sample_surface(950.0)   # deepest  (injection layer)
surf_L2 = sample_surface(900.0)
surf_L3 = sample_surface(850.0)

# We create a `GenericTopography` from the three surfaces, which implements the `AbstractTopography` interface required by `CO2BatchFill`. 
sand_layers = [
    Dict{String,Any}("name" => "Storage layer 3", "top" => surf_L3, "base" => surf_L3 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 2", "top" => surf_L2, "base" => surf_L2 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 1", "top" => surf_L1, "base" => surf_L1 .+ LAYER_THICK),
]
topo = GenericTopography(sand_layers, NX, NY, DX, DY, minimum(surf_L3), maximum(surf_L1) + LAYER_THICK)
domain = create_domain(topo, 0.1)

# Using the `spillanalysis` algorithm from SWIM, the `analyze_base_surfaces` function identifies the spill points and leakage paths for each layer. The `boundary_condition` argument controls whether $\text{CO}_2$ can leave the domain (open BC) or if boundary walls are added to prevent outflow (closed BC). 
layers = analyze_base_surfaces(topo; boundary_condition=boundary_condition)

# ## Defining injection schedule
TOTAL_RATE = 80_000.0
INJECTION_END = 10.0 # injection up to 10 years, then stop

null_rate = create_no_injection(layers)

injection_loc_1 = (div(NX, 2), div(NY, 2))
injection_loc_2 = (div(NX, 4), div(NY, 4))
rate_L1_1 = create_injection_rate(layers, injection_loc_1, TOTAL_RATE; radius=0)
rate_L1_2 = create_injection_rate(layers, injection_loc_2, TOTAL_RATE; radius=0)
rate_L1 = rate_L1_1 .+ rate_L1_2


# We create an `injection_events` vector with one entry per layer, where each entry is a vector of `InjectionEvent`s.
bottom_injection_events = [
    InjectionEvent(0.0, rate_L1),
    InjectionEvent(INJECTION_END, create_injection_rate(layers, (1, 1), 0.0)),
]


injection_events = [
    bottom_injection_events, # Layer 1
    null_rate, # Layer 2
    null_rate, # Layer 3
]

# ## Run simulation
# The `fill_layers` function runs the simulation. 
seqs, leakage_states, weather_events_per_layer = fill_layers(
    layers, domain, all_layers_rp, injection_events; verbose=false)

# After the simulation is complete, we can generate snapshots at arbitrary time points. `generate_multi_layer_snapshots` creates a vector of `MultiLayerSnapshot`s containing information about the total $\text{CO}_2$ distribution across all layers, and a `LayerSnapshot`s for each individual layer.
t_end = 15.0
timepoints = collect(range(0.0, stop=t_end, length=30))
multi_snaps = generate_multi_layer_snapshots(
    layers, seqs, leakage_states, weather_events_per_layer, timepoints)

# Print summary at final snapshot
print_summary(multi_snaps[end])

# using CairoMakie
# plot_multi_layer(
#     layers, multi_snaps[end], domain;
#     max_co2_height=ceil(round(rp.leakage_height, digits=2)),
#     show_contours=true,
#     show_labels=true,
#     contour_levels=20,
#     major_contour_every=5,
#     contour_opacity=1.0,
#     figure_size=(500 * N_LAYERS, 500),
#     colormap=:Blues,
#     injection_locations=[(Float64(injection_loc_1[1]) * DX, Float64(injection_loc_1[2]) * DY), (Float64(injection_loc_2[1]) * DX, Float64(injection_loc_2[2]) * DY)],
#     show_leakage_locations=true,
#     show_extents=false,
#     colorbar_label="Column height",
# )