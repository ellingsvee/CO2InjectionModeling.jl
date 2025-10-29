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
data = NPZ.npzread(datapath)
caprock_topography = data["topographies"][:, :, :]

# Load grid metadata for physical dimensions
xmin, xmax = data["xmin"], data["xmax"]
ymin, ymax = data["ymin"], data["ymax"]
length_x = Float64(xmax - xmin)  # Physical length in x-direction (meters)
length_y = Float64(ymax - ymin)  # Physical length in y-direction (meters)

println("Grid dimensions:")
println("  X range: [$xmin, $xmax] m, length = $(length_x) m")
println("  Y range: [$ymin, $ymax] m, length = $(length_y) m")

# ============================================================================
# 2. Analyze Layer Trap Structures
# ============================================================================

tstructs = analyze_layers(caprock_topography; lengths=(length_x, length_y))

# Select injection location
top_tstruct = tstructs[1]
# Select the injection location as the lowest point in the top layer
min_elevation = Inf
min_point = nothing
for footprint in top_tstruct.footprints
    elevations = top_tstruct.topography[footprint]
    # Find the point with the minimum elevation in this footprint
    local_min_elevation = minimum(elevations)
    if local_min_elevation < min_elevation
        min_elevation = local_min_elevation
        local_min_index = findfirst(==(local_min_elevation), elevations)
        min_point = footprint[local_min_index]
    end
end

injection_location = min_point

# Injection rate
nx, ny = size(caprock_topography[1, :, :])
ratio = nx * ny / (length_x * length_y)
injection_rate = 20000.0

# ============================================================================
# 2. Generate Ensemble of Leakage Heights
# ============================================================================

n_simulations = 25  # Number of ensemble members (increased for better statistics)

# # Generate different leakage heights by varying the fraction parameter
# fractions = range(1.0, 2.0, length=n_simulations)  # Different leakage thresholds

# leakage_height_samples = [
#     compute_leakage_heights(tstructs; fraction=frac)
#     for frac in fractions
# ]

heights = [40.8]
leakage_height_samples = [
    fill(height, length(tstructs) - 1) for height in heights
]


# ============================================================================
# 3. Run Ensemble Simulations
# ============================================================================

common_times, all_texs, all_n_filled, simulation_metadata = run_ensemble_simulations(
    tstructs,
    injection_location,
    injection_rate,
    leakage_height_samples;
    n_common_times=50,
    start_time=0.001
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