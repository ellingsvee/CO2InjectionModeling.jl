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

# Reservoir properties
rp = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, 25_000.0, 5.0)
rp_caprock = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, Inf, 5.0)

println("Leakage height threshold: $(round(rp.leakage_height, digits=2)) m")

# Injection — single central well in layer 1 (deepest) only
TOTAL_RATE = 80_000.0
INJECTION_END = 10.0

pad = PAD_WIDTH # since closed boundaries
topo_size = (NX + 2 * pad, NY + 2 * pad)

injection_location_idx = (div(NX, 2) + pad, div(NY, 2) + pad)
rate_L1 = zeros(topo_size)
rate_L1[injection_location_idx[1], injection_location_idx[2]] = TOTAL_RATE

injection_events = [
    # Layer 1: inject then stop
    [InjectionEvent(0.0, rate_L1), InjectionEvent(INJECTION_END, zeros(topo_size))],
    # Layers 2–3: no direct injection (CO2 arrives via leakage from below)
    [InjectionEvent(0.0, zeros(topo_size))],
    [InjectionEvent(0.0, zeros(topo_size))],
]

# Run multi-layer simulation
println("\nRunning multi-layer fill simulation...")
seqs, leakage_states, weather_events_per_layer = fill_layers(
    layers, domain, [rp, rp, rp_caprock], injection_events; verbose=false)

# Determine time range for snapshots
t_end = 15.0
timepoints = collect(range(0.0, stop=t_end, length=30))

println("\nGenerating $(length(timepoints)) snapshots from t=0 to t=$(round(t_end, digits=2))...")
multi_snaps = generate_multi_layer_snapshots(
    layers, seqs, leakage_states, weather_events_per_layer, timepoints)

# Print summary at final snapshot
println()
print_summary(multi_snaps[end])

# Settings for plotting
update_theme!(
    merge(
        theme_latexfonts(),
        Theme(fontsize=25)
    )
)

# Plot 1: static spatial distribution at final time (one panel per layer)
injection_location_loc = (div(NX, 2) * DX, div(NY, 2) * DY)
plot_multi_layer(
    layers, multi_snaps[end], domain;
    output_file=joinpath(@__DIR__, "multi_layer_co2_final.svg"),
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
)

# Plot 2: volume time-series per layer
_phys_scale = swim_volume_to_physical_volume(1.0, rp, domain)
plot_multi_layer_volumes_timeseries(
    multi_snaps;
    output_file=joinpath(@__DIR__, "multi_layer_timeseries_per_layer.svg"),
    vol_scale=_phys_scale / 1e5,
    ylabel=L"Volume $\left(\!\times\! 10^5\right)$",
)
