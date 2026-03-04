using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using CairoMakie
using Random

Random.seed!(42)

# ---------------------------------------------------------------------------
# Grid / domain settings
# ---------------------------------------------------------------------------
const NX, NY = 100, 100
const LENGTH_X = 1000.0   # m
const LENGTH_Y = 1000.0   # m
const DX, DY = LENGTH_X / NX, LENGTH_Y / NY
const LAYER_THICK = 10.0     # m (vertical thickness of each sand layer)
const PAD_WIDTH = 2
const N_LAYERS = 4

# ---------------------------------------------------------------------------
# Generate four GRF surfaces at increasing depth
# ---------------------------------------------------------------------------
cov = CovarianceFunction(2, Matern(75, 2, σ=1.0))
pts = range(0.0, stop=(NX - 1) * DX, step=DX)

function sample_surface(base_depth)
    sample(GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)) .+ base_depth
end

surf_L1 = sample_surface(950.0)   # deepest  (injection layer)
surf_L2 = sample_surface(900.0)
surf_L3 = sample_surface(850.0)
surf_L4 = sample_surface(800.0)   # shallowest

# ---------------------------------------------------------------------------
# Build topography — sand_layers ordered shallowest-first, deepest-last
# ---------------------------------------------------------------------------
sand_layers = [
    Dict{String,Any}("name" => "L4", "top" => surf_L4, "base" => surf_L4 .+ LAYER_THICK),
    Dict{String,Any}("name" => "L3", "top" => surf_L3, "base" => surf_L3 .+ LAYER_THICK),
    Dict{String,Any}("name" => "L2", "top" => surf_L2, "base" => surf_L2 .+ LAYER_THICK),
    Dict{String,Any}("name" => "L1", "top" => surf_L1, "base" => surf_L1 .+ LAYER_THICK),
]
topo = GenericTopography(sand_layers, NX, NY, DX, DY, minimum(surf_L4), maximum(surf_L1) + LAYER_THICK)
domain = create_domain(topo, 1.0)

# analyze_base_surfaces reverses the order → layers[1] = L1 (deepest, injection)
boundary_condition = :closed
layers = analyze_base_surfaces(topo; boundary_condition, pad_width=PAD_WIDTH)

# ---------------------------------------------------------------------------
# Reservoir properties
# ---------------------------------------------------------------------------
rp = ReservoirProperties(0.3, 0.2, 0.1, 10_000.0, 10.0)
rp_caprock = ReservoirProperties(0.3, 0.2, 0.1, Inf, 10.0)

println("Leakage height threshold: $(round(rp.leakage_height, digits=2)) m")

# ---------------------------------------------------------------------------
# Injection — single central well in layer 1 (deepest) only
# ---------------------------------------------------------------------------
TOTAL_RATE = 30_000.0
INJECTION_END = 10.0

pad = PAD_WIDTH # since closed boundaries
topo_size = (NX + 2 * pad, NY + 2 * pad)

rate_L1 = zeros(topo_size)
rate_L1[div(NX, 2)+pad, div(NY, 2)+pad] = TOTAL_RATE
# rate_L1[1+pad, 1+pad] = TOTAL_RATE / 2
# rate_L1[NX+pad, NY+pad] = TOTAL_RATE / 2

injection_events = [
    # Layer 1: inject then stop
    [InjectionEvent(0.0, rate_L1), InjectionEvent(INJECTION_END, zeros(topo_size))],
    # Layers 2–4: no direct injection (CO2 arrives via leakage from below)
    [InjectionEvent(0.0, zeros(topo_size))],
    [InjectionEvent(0.0, zeros(topo_size))],
    [InjectionEvent(0.0, zeros(topo_size))],
]

# ---------------------------------------------------------------------------
# Run multi-layer simulation
# ---------------------------------------------------------------------------
println("\nRunning multi-layer fill simulation...")
seqs, leakage_states, weather_events_per_layer = fill_layers(
    layers, domain, [rp, rp, rp, rp_caprock], injection_events; verbose=false)

# ---------------------------------------------------------------------------
# Determine time range for snapshots
# ---------------------------------------------------------------------------
timepoints = collect(range(0.0, stop=15.0, length=30))

println("\nGenerating $(length(timepoints)) snapshots from t=0 to t=$(round(t_end, digits=2))...")
multi_snaps = generate_multi_layer_snapshots(
    layers, seqs, leakage_states, weather_events_per_layer, timepoints)

# ---------------------------------------------------------------------------
# Print summary at final snapshot
# ---------------------------------------------------------------------------
println()
print_summary(multi_snaps[end])

# ---------------------------------------------------------------------------
# Plot 1: static spatial distribution at final time (one panel per layer)
# ---------------------------------------------------------------------------
println("\nPlotting spatial CO2 distribution...")
plot_multi_layer(
    layers, multi_snaps[end], domain;
    output_file=joinpath(@__DIR__, "multi_layer_co2_final.svg"),
    max_co2_height=5.0,
    show_contours=true,
    contour_levels=12,
)

# ---------------------------------------------------------------------------
# Plot 2: volume time-series per layer
# ---------------------------------------------------------------------------
println("Plotting per-layer volume time-series...")
plot_multi_layer_volumes_timeseries(
    multi_snaps;
    output_file=joinpath(@__DIR__, "multi_layer_timeseries_per_layer.svg"),
    mode=:per_layer,
)

# ---------------------------------------------------------------------------
# Plot 3: system-level time-series
# ---------------------------------------------------------------------------
println("Plotting system-level volume time-series...")
plot_multi_layer_volumes_timeseries(
    multi_snaps;
    output_file=joinpath(@__DIR__, "multi_layer_timeseries_system.svg"),
    mode=:system,
)

# ---------------------------------------------------------------------------
# Animate: N side-by-side panels
# ---------------------------------------------------------------------------
println("Animating multi-layer filling...")
animate_multi_layer_filling(
    layers, seqs, domain;
    output_file=joinpath(@__DIR__, "multi_layer_filling.gif"),
    num_frames=60,
    end_time=t_end,
    fps=15,
    max_co2_height=5.0,
    show_contours=true,
    contour_levels=12,
)

