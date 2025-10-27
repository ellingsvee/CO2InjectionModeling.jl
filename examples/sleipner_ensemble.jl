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
injection_rate = 20000.0

# ============================================================================
# 2. Generate Ensemble of Leakage Heights
# ============================================================================

n_simulations = 25  # Number of ensemble members (increased for better statistics)

# Generate different leakage heights by varying the fraction parameter
fractions = range(1.0, 2.0, length=n_simulations)  # Different leakage thresholds

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

fill_fractions, fill_uncertainties = compute_ensemble_statistics(all_texs)

# ============================================================================
# 5. Visualize Ensemble Probability Animation
# ============================================================================

println("\nCreating ensemble probability animation...")
fig_prob = create_ensemble_probability_animation(
    caprock_topography,
    fill_fractions,
    common_times,
    "../media/ensemble_probability.gif";
    title="Ensemble CO2 Probability",
    n_cols=3,
    framerate=5,
    co2_colormap=:viridis,
    figure_size=(1600, 1300),
    show_terrain=true
)

# ============================================================================
# 6. Visualize Individual Simulations
# ============================================================================

# # Optionally create animations for individual realizations
# for i in [1, div(n_simulations, 2), n_simulations]  # First, middle, last
#     println("\nCreating animation for simulation $i (fraction=$(round(fractions[i], digits=2)))...")
#     fig = create_layer_animation(
#         caprock_topography,
#         all_texs[i],
#         common_times,
#         "../media/ensemble_sim_$i.gif";
#         title="CO2 Migration - Realization $i (fraction=$(round(fractions[i], digits=2)))",
#         n_cols=3,
#         framerate=5,
#         co2_colormap=:hot,
#         figure_size=(1600, 1300),
#         show_terrain=true
#     )
# end
# ============================================================================