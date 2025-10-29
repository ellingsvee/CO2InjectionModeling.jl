# ============================================================================
# Sleipner Ensemble Simulations with Extended Domain Boundary
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
caprock_topography_original = data["topographies"][:, :, :]

# Load grid metadata for physical dimensions
xmin, xmax = data["xmin"], data["xmax"]
ymin, ymax = data["ymin"], data["ymax"]

# Calculate physical dimensions
length_x = xmax - xmin
length_y = ymax - ymin

# ============================================================================
# 2. Extend Domain with Dome Boundary
# ============================================================================

# Add a boundary region that gradually rises to form a dome
# This prevents CO2 from flowing out of the domain unrealistically
boundary_width = 100  # Number of cells to add on each side (adjust as needed)

caprock_topography, inner_domain = extend_domain_with_dome(
    caprock_topography_original, 
    boundary_width;
    transition_type=:exponential,  # :linear, :exponential, or :sigmoid
    steepness=1.0  # Controls how quickly boundary rises (higher = steeper)
)

println("Original domain size: ", size(caprock_topography_original))
println("Extended domain size: ", size(caprock_topography))
println("Inner domain indices: ", inner_domain)

# Update physical dimensions to account for extended domain
# Assuming same cell size
nx_orig, ny_orig = size(caprock_topography_original[1, :, :])
cell_size_x = length_x / ny_orig
cell_size_y = length_y / nx_orig

# New physical dimensions
length_x_extended = length_x + 2 * boundary_width * cell_size_x
length_y_extended = length_y + 2 * boundary_width * cell_size_y

# ============================================================================
# 2.1 Visualize Extended Domain
# ============================================================================

println("\nCreating domain extension visualization...")

# Select a representative layer (e.g., first layer) to visualize
layer_idx = 2

fig_domain = Figure(size=(1400, 600))

# Original domain
ax1 = Axis(fig_domain[1, 1];
    aspect=DataAspect(),
    title="Original Domain - Layer $layer_idx",
    xlabel="Y index",
    ylabel="X index")

original_layer = reverse(transpose(caprock_topography_original[layer_idx, :, :]), dims=1)
hm1 = heatmap!(ax1, original_layer; colormap=:terrain)
Colorbar(fig_domain[1, 2], hm1, label="Elevation (m)")

# Extended domain with boundary highlighted
ax2 = Axis(fig_domain[1, 3];
    aspect=DataAspect(),
    title="Extended Domain - Layer $layer_idx (boundary highlighted)",
    xlabel="Y index",
    ylabel="X index")

extended_layer = reverse(transpose(caprock_topography[layer_idx, :, :]), dims=1)
hm2 = heatmap!(ax2, extended_layer; colormap=:terrain)
Colorbar(fig_domain[1, 4], hm2, label="Elevation (m)")

# Draw rectangle showing original domain boundary
# Recompute visualization coords with correct axis mapping (x <- j, y <- i)
ny_vis, nx_vis = size(extended_layer)
i_range, j_range = inner_domain

vis_x_start = minimum(j_range)
vis_x_end   = maximum(j_range)
vis_y_start = minimum(i_range) 
vis_y_end   = maximum(i_range)

lines!(ax2, [vis_x_start, vis_x_end, vis_x_end, vis_x_start, vis_x_start],
           [vis_y_start, vis_y_start, vis_y_end, vis_y_end, vis_y_start];
           color=:red, linewidth=3, label="Original domain")
axislegend(ax2, position=:rt)

save("../media/domain_extension_visualization.png", fig_domain)
println("Domain visualization saved to ../media/domain_extension_visualization.png")
display(fig_domain)


# ============================================================================
# 3. Analyze Layer Trap Structures (on extended domain)
# ============================================================================

tstructs = analyze_layers(caprock_topography; lengths=(length_x_extended, length_y_extended))

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
ratio = nx * ny / (length_x_extended * length_y_extended)
injection_rate = 20000.0

# ============================================================================
# 4. Generate Ensemble of Leakage Heights
# ============================================================================

n_simulations = 25  # Number of ensemble members

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
# 5. Run Ensemble Simulations
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
# 6. Compute Ensemble Statistics
# ============================================================================

fill_fractions, fill_uncertainties = compute_ensemble_statistics(all_texs)

# ============================================================================
# 7. Visualize Ensemble Probability Animation (showing only inner domain)
# ============================================================================

println("\nCreating ensemble probability animation (inner domain only)...")
fig_prob = create_ensemble_probability_animation(
    caprock_topography,
    fill_fractions,
    common_times,
    "../media/ensemble_probability_extended_domain.gif";
    title="Ensemble CO2 Probability (Extended Domain)",
    n_cols=3,
    framerate=5,
    co2_colormap=:viridis,
    figure_size=(1600, 1300),
    show_terrain=true,
    inner_domain=inner_domain  # Only visualize the inner (original) domain
)

println("Simulation complete!")
println("Extended domain prevents CO2 from flowing out while maintaining natural geometry.")
println("Visualization shows only the original (inner) domain.")
