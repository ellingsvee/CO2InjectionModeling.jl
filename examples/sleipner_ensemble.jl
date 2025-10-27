# ============================================================================
# Sleipner Ensemble Simulations - Uncertainty Quantification
# ============================================================================

using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using CairoMakie
using ColorSchemes
using NPZ
using Random
using Statistics

# Import RGBA for probability visualization
using Colors: RGBA

# ============================================================================
# 1. Load Data and Setup
# ============================================================================

datapath = "../sleipner/sleipner_topographies.npz"
caprock_topography = NPZ.npzread(datapath)["topographies"][:, :, :]

tstructs = analyze_layers(caprock_topography)

# Select injection location
top_tstruct = tstructs[1]
footprints = top_tstruct.footprints
trap_idx = rand(1:length(footprints))
footprint = footprints[trap_idx]
injection_location = rand(footprint)
injection_rate = 1.0

# ============================================================================
# 2. Generate Ensemble of Leakage Heights
# ============================================================================

n_simulations = 10  # Number of ensemble members (increased for better statistics)

# Generate different leakage heights by varying the fraction parameter
fractions = range(1.2, 2.0, length=n_simulations)  # Different leakage thresholds

leakage_height_samples = [
    compute_leakage_heights(tstructs; fraction=frac)
    for frac in fractions
]

println("Running $n_simulations ensemble simulations with fraction range: $(fractions[1]) to $(fractions[end])")

# ============================================================================
# 3. Run Ensemble Simulations
# ============================================================================

common_times, all_texs, all_n_filled, simulation_metadata = run_ensemble_simulations(
    tstructs,
    injection_location,
    injection_rate,
    leakage_height_samples;
    n_common_times=100,
    start_time=0.1
)

# ============================================================================
# 4. Compute Ensemble Statistics
# ============================================================================

println("\nComputing ensemble statistics...")

# Additional diagnostic: Check structure of all_texs
println("\nDiagnostic - all_texs structure:")
println("  Number of simulations: $(length(all_texs))")
for sim_idx in 1:length(all_texs)
    println("  Simulation $sim_idx: $(length(all_texs[sim_idx])) layers")
    for layer_idx in 1:min(3, length(all_texs[sim_idx]))
        println("    Layer $layer_idx: $(length(all_texs[sim_idx][layer_idx])) time steps")
        # Check a sample state value
        sample_state = all_texs[sim_idx][layer_idx][1]
        n_co2 = count(x -> x > 1, sample_state)
        println("      Time 1: CO2 cells = $n_co2 / $(length(sample_state))")
    end
end

fill_fractions, fill_uncertainties = compute_ensemble_statistics(all_texs)

# Diagnostic: Check if we have non-zero probabilities
println("\nDiagnostic - Checking fill_fractions data:")
for layer_idx in 1:min(3, length(fill_fractions))
    for time_idx in [1, div(length(fill_fractions[layer_idx]), 2), length(fill_fractions[layer_idx])]
        max_prob = maximum(fill_fractions[layer_idx][time_idx])
        min_prob = minimum(fill_fractions[layer_idx][time_idx])
        n_nonzero = count(x -> x > 0, fill_fractions[layer_idx][time_idx])
        println("  Layer $layer_idx, Time $time_idx: min=$(round(min_prob, digits=3)), max=$(round(max_prob, digits=3)), nonzero=$n_nonzero")
    end
end

println("\nEnsemble Summary:")
for (i, meta) in enumerate(simulation_metadata)
    println("Simulation $i: reached $(meta.n_filled_layers) layers, fraction=$(round(fractions[i], digits=2))")
end

# ============================================================================
# 5. Visualize Ensemble Probability Animation
# ============================================================================

println("\nCreating ensemble probability animation...")
fig_prob = create_ensemble_probability_animation(
    caprock_topography,
    fill_fractions,
    common_times,
    "../media/ensemble_probability.gif";
    title="Ensemble CO2 Probability (Opacity = Fill Fraction)",
    n_cols=3,
    framerate=5,
    topo_colormap=:terrain,
    co2_colormap=:hot,
    figure_size=(1600, 1300),
    show_terrain=false,
    terrain_contours=false
)

# ============================================================================
# 6. Visualize Individual Simulations
# ============================================================================

# Optionally create animations for individual realizations
for i in [1, 3, n_simulations]  # First, middle, last
    println("\nCreating animation for simulation $i (fraction=$(round(fractions[i], digits=2)))...")
    fig = create_layer_animation(
        caprock_topography,
        all_texs[i],
        common_times,
        "../media/ensemble_sim_$i.gif";
        title="CO2 Migration - Realization $i (fraction=$(round(fractions[i], digits=2)))",
        n_cols=3,
        framerate=5,
        topo_colormap=:terrain,
        co2_colormap=:hot,
        figure_size=(1600, 1300)
    )
end

# ============================================================================
# 7. Analyze Uncertainty at Specific Time
# ============================================================================

# Choose a specific time point (e.g., halfway through simulation)
time_idx = div(length(common_times), 2)
time_point = common_times[time_idx]

println("\nAnalyzing uncertainty at t=$(round(time_point, digits=2)) years...")

# Create figure showing fill probability and uncertainty for each layer at this time
fig_uncertainty = Figure(size=(1800, 1200))

n_layers = length(tstructs)
n_cols = 3
n_rows = ceil(Int, n_layers / n_cols)

for layer_idx in 1:n_layers
    row = div(layer_idx - 1, n_cols) + 1
    col = mod(layer_idx - 1, n_cols) + 1

    # Fill probability (fraction of simulations with CO2)
    ax_prob = Axis(fig_uncertainty[2*row-1, col];
        aspect=DataAspect(),
        title="Layer $layer_idx - Fill Probability",
        titlesize=14)

    fill_prob = reverse(transpose(fill_fractions[layer_idx][time_idx]), dims=1)
    heatmap!(ax_prob, fill_prob; colormap=:viridis, colorrange=(0, 1))
    hidedecorations!(ax_prob)

    # Uncertainty (standard deviation of fill indicator)
    ax_std = Axis(fig_uncertainty[2*row, col];
        aspect=DataAspect(),
        title="Layer $layer_idx - Uncertainty",
        titlesize=14)

    uncertainty = reverse(transpose(fill_uncertainties[layer_idx][time_idx]), dims=1)
    heatmap!(ax_std, uncertainty; colormap=:reds, colorrange=(0, 0.5))
    hidedecorations!(ax_std)
end

Label(fig_uncertainty[0, :], "Ensemble Statistics at t=$(round(time_point, digits=2)) years",
    fontsize=20, halign=:center, font=:bold)

save("../media/ensemble_uncertainty.png", fig_uncertainty)
println("Uncertainty plot saved to ../media/ensemble_uncertainty.png")

println("\nEnsemble analysis complete!")
