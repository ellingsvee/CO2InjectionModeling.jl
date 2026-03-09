using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using CairoMakie
using LaTeXStrings
using Random

Random.seed!(101)

# Grid / domain settings
const NX, NY = 100, 100
const LENGTH_X = 1000.0
const LENGTH_Y = 1000.0
const DX, DY = LENGTH_X / NX, LENGTH_Y / NY
const LAYER_THICK = 10.0     # m (vertical thickness of each sand layer)
const PAD_WIDTH = 2
const N_LAYERS = 3
const RESIDUAL_TRAPPING = 0.4

# Generate four GRF surfaces at increasing depth
cov = CovarianceFunction(2, Matern(200, 2, σ=3.0))
pts = range(0.0, stop=(NX - 1) * DX, step=DX)

function sample_surface(base_depth)
    sample(GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)) .+ base_depth
end

surf_L1 = sample_surface(950.0)   # deepest  (injection layer)
surf_L2 = sample_surface(900.0)
surf_L3 = sample_surface(850.0)

# Build topography — sand_layers ordered shallowest-first, deepest-last
sand_layers = [
    Dict{String,Any}("name" => "Storage layer 3", "top" => surf_L3, "base" => surf_L3 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 2", "top" => surf_L2, "base" => surf_L2 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 1", "top" => surf_L1, "base" => surf_L1 .+ LAYER_THICK),
]
topo = GenericTopography(sand_layers, NX, NY, DX, DY, minimum(surf_L3), maximum(surf_L1) + LAYER_THICK)
domain = create_domain(topo, 0.1)

# analyze_base_surfaces reverses the order → layers[1] = L1 (deepest, injection)
boundary_condition = :closed
layers = analyze_base_surfaces(topo; boundary_condition, pad_width=PAD_WIDTH)

# Injection — single central well in layer 1 (deepest) only
TOTAL_RATE = 80_000.0
INJECTION_END = 10.0

# create_injection_rate handles padding offset automatically
rate_L1 = create_injection_rate(layers, (div(NX, 2), div(NY, 2)), TOTAL_RATE)

injection_events = [
    # Layer 1: inject then stop
    [InjectionEvent(0.0, rate_L1), InjectionEvent(INJECTION_END, create_injection_rate(layers, (1, 1), 0.0))],
    # Layers 2–3: no direct injection (CO2 arrives via leakage from below)
    create_no_injection(layers),
    create_no_injection(layers),
]

rp_caprock = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, Inf, 5.0)

N_ENSEMBLE = 100
capillary_entry_pressures = range(20_000.0, stop=30_000.0, length=N_ENSEMBLE)
println("Min meakage height threshold: $(round(ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, minimum(capillary_entry_pressures), 5.0).leakage_height, digits=2))")
println("Max meakage height threshold: $(round(ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, maximum(capillary_entry_pressures), 5.0).leakage_height, digits=2))")

t_end = 15.0
timepoints = collect(range(0.0, stop=t_end, length=30))

multi_snaps_ensemble_final = Vector{MultiLayerSnapshot}(undef, N_ENSEMBLE)
multi_snaps_ensemble_ts = Vector{Vector{MultiLayerSnapshot}}(undef, N_ENSEMBLE)

# Time the ensemble runs
@time for i in 1:N_ENSEMBLE
    rp = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, capillary_entry_pressures[i], 5.0)

    # Run multi-layer simulation
    seqs, leakage_states, weather_events_per_layer = fill_layers(
        layers, domain, [rp, rp, rp_caprock], injection_events; verbose=false)

    multi_snaps = generate_multi_layer_snapshots(
        layers, seqs, leakage_states, weather_events_per_layer, timepoints)

    multi_snaps_ensemble_final[i] = multi_snaps[end]
    multi_snaps_ensemble_ts[i] = multi_snaps
end

# Settings for plotting
update_theme!(
    merge(
        theme_latexfonts(),
        Theme(fontsize=25)
    )
)

# Plot ensemble volume timeseries — mean ± 1σ bands per layer
_phys_scale = volume_scale(ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, 1.0, 5.0), domain)
plot_multi_layer_ensemble_timeseries(
    multi_snaps_ensemble_ts;
    output_file=joinpath(@__DIR__, "ensemble_timeseries.svg"),
    vol_scale=_phys_scale / 1e5,
    ylabel=L"Volume $\left(\!\times\! 10^5\right)$",
)
