
using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using CairoMakie
using Printf

include("monte_carlo_analysis.jl")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Create Monte Carlo configuration
# You can customize parameter ranges here
# config = create_custom_monte_carlo_config(
#     # Number of Monte Carlo realizations
#     n_realizations=10,

#     # Parameter ranges (min, max) for each reservoir property
#     sand_porosity=(0.32, 0.48),                          # ±20% around 0.4
#     sand_residual_co2_saturation=(0.10, 0.40),           # ±25% around 0.2
#     sand_irreducible_water_saturation=(0.20, 0.34),      # ±20% around 0.3
#     shale_pressure_threshold=(80000.0, 116000.0),        # ±18% around 98000 Pa
#     # leakage_height=(16.0, 25.0),                         # ±20% around 20 m
#     # residual_leakage_time=(3.0, 7.0),                    # ±40% around 5 years
#     residual_leakage_time=(2.0, 8.0),                    # ±40% around 5 years

#     # Simulation time settings
#     start_time=0.0,
#     end_time=15.0,
#     num_snapshots=30,

#     # Random seed for reproducibility (set to nothing for random)
#     random_seed=42
# )

 config = create_default_monte_carlo_config(
    n_realizations=25,
    start_time=0.0,
    end_time=15.0,
    num_snapshots=45,
    random_seed=1
)

# =============================================================================
# SETUP RESERVOIR GEOMETRY AND INJECTION
# =============================================================================

# Load Sleipner topography and create domain
boundary_condition = :closed
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)


# plot_layer_topographies(
#     layers,
#     domain;
#     output_file="layer_topographies.png",
#     title="Sleipner Layer Topographies",
#     colormap=:viridis,
#     contour_levels=14,
#     show_contour_labels=false
# );

# Load injection location from feeder chimney data
xy, (utm_x, utm_y, depth) = load_feeder_location(topography)


# Generate injection schedule
injection_events = generate_sleipner_injection_events(layers, xy)

# =============================================================================
# RUN MONTE CARLO SIMULATION
# =============================================================================

ensemble = run_monte_carlo_simulation(
    config,
    layers,
    domain,
    injection_events;
    verbose=true,
    store_seqs=true  # Store sequences for spatial animations
);


# =============================================================================
# UNCERTAINTY ANALYSIS
# =============================================================================

# Create output directory for plots
mkpath("monte_carlo_results")

# # 1. Total CO2 storage uncertainty
# plot_total_storage_uncertainty(
#     ensemble;
#     output_file="monte_carlo_results/total_storage_uncertainty.png",
#     show_individual=false,
#     figsize=(1000, 600)
# )

# # 2. Total CO2 leakage uncertainty
# println("\nPlotting total CO2 leakage uncertainty...")
# plot_total_leakage_uncertainty(
#     ensemble;
#     output_file="monte_carlo_results/total_leakage_uncertainty.png",
#     show_individual=false,
#     figsize=(1000, 600)
# )

# # 3. Layer-by-layer uncertainty
# println("\nPlotting layer-by-layer uncertainty...")
# n_layers = length(layers)
# plot_all_layers_uncertainty(
#     ensemble,
#     n_layers;
#     output_file="monte_carlo_results/all_layers_uncertainty.png",
#     figsize=(1400, 1000)
# )

# # 4. Layer distribution barplot at final timepoint
# println("\nGenerating layer distribution barplot...")
plot_layer_distribution_barplot(
    ensemble;
    output_file="monte_carlo_results/layer_distribution_barplot.svg",
    bar_color=colorant"#386624",
    bar_alpha=1.0,  # Semi-transparent bars so error bars are visible
    error_bar_color=:black,  # Black error bars for better contrast
    error_bar_linewidth=2.0,
    fontsize_labels=16,
    fontsize_ticks=15,
    fontsize_values=15,
    figure_size=(800, 400)
)
println("jk")

# # 5. Mean plume animation across realizations (single color)
# println("\nGenerating mean plume animation (single color)...")
# animate_multi_layer_filling_ensemble(
#     ensemble,
#     layers,
#     domain;
#     output_file="monte_carlo_results/ensemble_mean_filling.gif",
#     num_frames=30,
#     end_time=15.0,
#     fps=3,
#     show_probability=false,
#     plume_color=:red
# )

# 6. Probability-based plume animation across realizations
        # :shale_pressure_threshold_std => Uniform(8.0, 15.0)
# animate_multi_layer_filling_ensemble(
#     ensemble,
#     layers,
#     domain;
#     output_file="monte_carlo_results/ensemble_probability_filling.gif",
#     num_frames=30,
#     end_time=15.0,
#     fps=3,
#     show_probability=true,
#     probability_colormap=:viridis
# )

# =============================================================================
# STATISTICAL SUMMARY
# =============================================================================

# Print comprehensive statistical summary
# print_uncertainty_summary(ensemble)
print_percentages(ensemble);

