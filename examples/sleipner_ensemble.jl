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
# boundary_width = 50
boundary_width = 5

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
area_extended = length_x_extended * length_y_extended

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

injection_rates_per_year_Mt = [0.07, 0.67, 0.85, 0.94, 0.94, 1.02, 0.96, 0.92, 0.76, 0.87, 0.83, 0.93, 0.82, 0.86, 0.76] # From 1996 to 2010
# injection_rate_per_year_Mt = mean(injection_rates_per_year_Mt)
injection_rate_per_year_Mt = 1.0
# Convert from Mt/year to m³/year (assuming CO2 density ~650 kg/m³)
co2_density = 650.0  # kg/m³
injection_rate = (injection_rate_per_year_Mt .* 1e9) ./ co2_density  # Mt to kg, then to m³

# ============================================================================
# 4. Generate Ensemble of Leakage Heights
# ============================================================================

n_simulations = 25

function sample_leakage_heights(
    n_barriers::Int,  # Number of barriers (8 for 9 layers)
    n_samples::Int;
    capillary_threshold_means::Vector{Float64} = fill(0.05, n_barriers),  # Layer-specific means in MPa
    capillary_threshold_stds::Vector{Float64} = fill(0.016, n_barriers),  # Layer-specific stds in MPa
    brine_to_co2_difference::Float64 = 300.0,
    gravitational_constant::Float64 = 9.81,
)::Vector{Vector{Float64}}
    # Convert to Pa
    capillary_threshold_means_pa = capillary_threshold_means .* 1e6
    capillary_threshold_stds_pa = capillary_threshold_stds .* 1e6
    
    leakage_height_samples = Vector{Vector{Float64}}(undef, n_samples)
    
    for i in 1:n_samples
        sampled_thresholds = Vector{Float64}(undef, n_barriers)
        for j in 1:n_barriers
            # Sample from log-normal for realism (positive values, skewed)
            # LogNormal params: mean of log = log(mean), std of log = std/mean (approx CV)
            log_mean = log(capillary_threshold_means_pa[j])
            log_std = capillary_threshold_stds_pa[j] / capillary_threshold_means_pa[j]  # CV ~ std/mean
            sampled_thresholds[j] = rand(LogNormal(log_mean, log_std))
        end
        # Compute leakage heights
        leakage_heights = sampled_thresholds ./ (brine_to_co2_difference * gravitational_constant)
        leakage_height_samples[i] = leakage_heights
    end
    return leakage_height_samples
end

# Define layer-specific means: lower for thin shales (1-7), higher for thick (8)
n_barriers = length(tstructs) - 1  # Assuming 8
capillary_threshold_means_kPa = [100, 90, 95, 105, 90, 112, 75, 120]  # kPa
capillary_threshold_means = capillary_threshold_means_kPa ./ 1000.0  # Convert to MPa
capillary_threshold_stds = fill(0.016, n_barriers)  # MPa; uniform small variation

leakage_height_samples = sample_leakage_heights(
    n_barriers, n_simulations;
    capillary_threshold_means=capillary_threshold_means,
    capillary_threshold_stds=capillary_threshold_stds,
)

# ============================================================================
# 5. Run Ensemble Simulations
# ============================================================================

# Set maximum simulation time (in years) - comment out or set to nothing for full simulation
max_simulation_time = 30.0
# max_simulation_time = nothing  # Uncomment to run until last update

common_times, all_texs, all_n_filled, simulation_metadata = run_ensemble_simulations(
    tstructs,
    injection_location,
    injection_rate,
    leakage_height_samples;
    n_common_times=50,
    start_time=0.001,
    max_time=max_simulation_time  # Limit simulation to max_simulation_time years
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
    title="Ensemble CO2 Probability",
    n_cols=3,
    framerate=5,
    co2_colormap=:viridis,
    figure_size=(1600, 1300),
    show_terrain=true,
    inner_domain=inner_domain  # Only visualize the inner (original) domain
)

# ============================================================================
# 8. Visualize CO2 Volume Over Time
# ============================================================================

println("\nCreating CO2 volume timeseries plot...")

# Calculate cell volume (m³)
cell_volume = (length_x_extended / ny) * (length_y_extended / nx)

# Compute volumes for all simulations
all_volumes = compute_co2_volumes_over_time(tstructs, all_texs, common_times, cell_volume)

# Create volume timeseries plot with all layers
fig_vol = create_volume_timeseries_plot(
    tstructs,
    all_texs,
    common_times,
    cell_volume,
    "../media/co2_volume_timeseries.png";
    title="CO2 Storage Volume Over Time",
    figure_size=(1400, 1200),
    show_total=true,
    show_cumulative=true
)
