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
using Distributions

# Import RGBA for probability visualization
using Colors: RGBA

# ============================================================================
# 1. Load Data and Setup
# ============================================================================

datapath = "../sleipner/sleipner_topographies.npz"
data = NPZ.npzread(datapath)
caprock_topography_raw = data["topographies"][:, :, :]

# Load grid metadata for physical dimensions
xmin, xmax = data["xmin"], data["xmax"]
ymin, ymax = data["ymin"], data["ymax"]

# Add border padding to prevent CO2 from leaking out of domain
border_width = 1  # Number of cells to add on each side
n_layers, ny_original, nx_original = size(caprock_topography_raw)

# Create padded array with border
ny_padded = ny_original + 2 * border_width
nx_padded = nx_original + 2 * border_width
caprock_topography = zeros(Float64, n_layers, ny_padded, nx_padded)

# Set border to very large depth (much deeper than any real depth)
# This creates a "wall" that CO2 cannot cross
border_depth = maximum(caprock_topography_raw)

for layer_idx in 1:n_layers
    # Fill entire layer with border depth first
    caprock_topography[layer_idx, :, :] .= border_depth
    
    # Copy original data into the center region
    caprock_topography[layer_idx, 
                      (border_width+1):(ny_padded-border_width),
                      (border_width+1):(nx_padded-border_width)] = caprock_topography_raw[layer_idx, :, :]
end

# Update grid dimensions to account for padding
# Calculate cell sizes
dx_original = Float64(xmax - xmin) / nx_original
dy_original = Float64(ymax - ymin) / ny_original

# Extend the domain by the border width
xmin_padded = xmin - border_width * dx_original
xmax_padded = xmax + border_width * dx_original
ymin_padded = ymin - border_width * dy_original
ymax_padded = ymax + border_width * dy_original

length_x = Float64(xmax_padded - xmin_padded)
length_y = Float64(ymax_padded - ymin_padded)

println("Grid dimensions (with border padding):")
println("  Original size: ($ny_original, $nx_original)")
println("  Padded size: ($ny_padded, $nx_padded)")
println("  Border width: $border_width cells")
println("  X range: [$xmin_padded, $xmax_padded] m, length = $(length_x) m")
println("  Y range: [$ymin_padded, $ymax_padded] m, length = $(length_y) m")

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

function sample_leakage_heights(
    shape::Tuple,
    n_samples::Int;
    capillary_threshold_mean::Float64 = 0.1, # MPa
    capillary_threshold_std::Float64 = 0.03, # MPa
    brine_to_co2_difference::Float64 = 300.0,
    gravitational_constant::Float64 = 9.81,
)::Vector{Vector{Float64}}
    capillary_threshold_mean *= 1e6  # Convert MPa to Pa
    capillary_threshold_std *= 1e6   # Convert MPa to Pa

    leakage_height_samples = Vector{Vector{Float64}}(undef, n_samples)

    # Generate n_samples with the specified shape
    for i in 1:n_samples
        # Sample capillary thresholds from a normal distribution
        sampled_thresholds = rand(Normal(capillary_threshold_mean, capillary_threshold_std), shape)
        # Compute leakage heights
        leakage_heights = sampled_thresholds ./ (brine_to_co2_difference * gravitational_constant)
        leakage_height_samples[i] = vec(leakage_heights)
    end

    return leakage_height_samples
    
end

leakage_height_samples = sample_leakage_heights(
    (length(tstructs) - 1,), n_simulations;
    capillary_threshold_mean=0.1,
    capillary_threshold_std=0.03,
)


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

# Remove border padding for visualization (show only the interior domain)
caprock_topography_interior = caprock_topography[:, 
                                                  (border_width+1):(ny_padded-border_width),
                                                  (border_width+1):(nx_padded-border_width)]

# Remove padding from fill_fractions
# fill_fractions is a vector of length n_layers, each containing a vector of time steps,
# each time step containing a (ny_padded, nx_padded) array
fill_fractions_interior = [
    [timestep[(border_width+1):(ny_padded-border_width), 
              (border_width+1):(nx_padded-border_width)] 
     for timestep in layer_timeseries]
    for layer_timeseries in fill_fractions
]

println("\nCreating ensemble probability animation...")
fig_prob = create_ensemble_probability_animation(
    caprock_topography_interior,
    fill_fractions_interior,
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