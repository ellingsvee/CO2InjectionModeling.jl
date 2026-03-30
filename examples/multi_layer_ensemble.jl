# # Multi-layer ensemble
#
# This example demonstrates how the speed of `CO2BatchFill` enables large ensembles of multi-layer simulations. We use the same setup as in the multi-layer filling example. However, instead of a single run, we generate an ensemble of 100 simulations where the `CAPILLARY_ENTRY_PRESSURE` is varied.

# ## Setup 
# This is essentially the same as in the multi-layer filling example, so we won't go into as much detail here. 
using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields # for generating synthetic surfaces
using CairoMakie # for visualization
using LaTeXStrings
using Random

# Physical properties and domain setup
const N_LAYERS = 3
const NX, NY = 100, 100
const LENGTH_X = 1000.0
const LENGTH_Y = 1000.0
const DX, DY = LENGTH_X / NX, LENGTH_Y / NY
const LAYER_THICK = 10.0     # m (vertical thickness of each sand layer)
const RESIDUAL_TRAPPING = 0.4 # fraction of CO2 volume that becomes immobile when a layer fills
const boundary_condition = :closed

# Set up the properties for the caproch layer
rp_caprock = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, Inf, 5.0) # Infinite entry pressure for caprock layer 


# Generating synthetic topography
Random.seed!(101) # for reproducibility
cov = CovarianceFunction(2, Matern(200, 2, σ=3.0))
pts = range(0.0, stop=(NX - 1) * DX, step=DX)

function sample_surface(base_depth)
    sample(GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)) .+ base_depth
end

surf_L1 = sample_surface(950.0)   # deepest  (injection layer)
surf_L2 = sample_surface(900.0)
surf_L3 = sample_surface(850.0)

# Setting up and analyzing the layers
sand_layers = [
    Dict{String,Any}("name" => "Storage layer 3", "top" => surf_L3, "base" => surf_L3 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 2", "top" => surf_L2, "base" => surf_L2 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 1", "top" => surf_L1, "base" => surf_L1 .+ LAYER_THICK),
]
topo = GenericTopography(sand_layers, NX, NY, DX, DY, minimum(surf_L3), maximum(surf_L1) + LAYER_THICK)
domain = create_domain(topo, 0.1)
layers = analyze_base_surfaces(topo; boundary_condition=boundary_condition)

# Defining injection schedule
TOTAL_RATE = 80_000.0
INJECTION_END = 10.0 # injection up to 10 years, then stop
rate_L1 = create_injection_rate(layers, (div(NX, 2), div(NY, 2)), TOTAL_RATE)
null_rate = create_no_injection(layers)

# Creating an `injection_events` vector with one entry per layer, where each entry is a vector of `InjectionEvent`s.
injection_events = [
    [InjectionEvent(0.0, rate_L1), InjectionEvent(INJECTION_END, create_injection_rate(layers, (1, 1), 0.0))], # Layer 1
    null_rate, # Layer 2
    null_rate, # Layer 3
]

# ## Setting up and running ensemble
# We create an ensemble of 100 simulations where the `CAPILLARY_ENTRY_PRESSURE` is varied between $20$ and $30$ kPa
N_ENSEMBLE = 100
capillary_entry_pressures = range(20_000.0, stop=30_000.0, length=N_ENSEMBLE)
println("Min meakage height threshold: $(round(ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, minimum(capillary_entry_pressures), 5.0).leakage_height, digits=2))")
println("Max meakage height threshold: $(round(ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, maximum(capillary_entry_pressures), 5.0).leakage_height, digits=2))")


# Empty vectors to store snapshots for each ensemble member
t_end = 15.0
timepoints = collect(range(0.0, stop=t_end, length=30))
multi_snaps_ensemble_final = Vector{MultiLayerSnapshot}(undef, N_ENSEMBLE)
multi_snaps_ensemble_ts = Vector{Vector{MultiLayerSnapshot}}(undef, N_ENSEMBLE)

# The `@time` macro is used to measure the total time taken for the ensemble simulations. 
@time for i in 1:N_ENSEMBLE
    ## Creating temporary ReservoirProperties struct
    rp = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, capillary_entry_pressures[i], 5.0)

    ## Running simulation
    seqs, leakage_states, weather_events_per_layer = fill_layers(
        layers, domain, [rp, rp, rp_caprock], injection_events; verbose=false)

    ## Generate snapshots for this ensemble member
    multi_snaps = generate_multi_layer_snapshots(
        layers, seqs, leakage_states, weather_events_per_layer, timepoints)

    ## Store the final snapshot and the full timeseries of snapshots for this ensemble member
    multi_snaps_ensemble_final[i] = multi_snaps[end]
    multi_snaps_ensemble_ts[i] = multi_snaps
end

# ## Visualization
# Optional defaults for making the plots look nicer.
update_theme!(
    merge(
        theme_latexfonts(),
        Theme(fontsize=25)
    )
)


# The `plot_multi_layer_ensemble` function plots the probability of the plume extents over the ensemble
injection_location_loc = (div(NX, 2) * DX, div(NY, 2) * DY)
plot_multi_layer_ensemble(
    layers, multi_snaps_ensemble_final, domain;
    show_contours=true,
    show_labels=true,
    contour_levels=20,
    major_contour_every=5,
    contour_opacity=1.0,
    figure_size=(500 * N_LAYERS, 500),
    colormap=:viridis,
    injection_locations=[injection_location_loc],
    show_leakage_locations=true,
)

# The `plot_multi_layer_ensemble_timeseries` is similar to `plot_multi_layer_timeseries`, but it plots the means and standard deviations across the ensemble.
_phys_scale = volume_scale(ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, 1.0, 5.0), domain)
plot_multi_layer_ensemble_timeseries(
    multi_snaps_ensemble_ts;
    output_file=joinpath(@__DIR__, "ensemble_timeseries.svg"),
    vol_scale=_phys_scale / 1e5,
    ylabel=L"Volume $\left(\!\times\! 10^5\right)$",
)
