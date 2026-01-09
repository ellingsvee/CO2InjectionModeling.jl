"""
Monte Carlo Uncertainty Quantification for Sleipner CO2 Injection

This script performs Monte Carlo simulations to quantify uncertainty in CO2 storage
predictions arising from uncertain reservoir properties.

The analysis produces:
1. Uncertainty bands for total CO2 storage and leakage over time
2. Layer-by-layer uncertainty analysis
3. Statistical summaries at key time points
4. Comprehensive visualizations

Usage:
    julia --project=. scripts/monte_carlo_simulation.jl
"""


using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using Printf

include("monte_carlo_analysis.jl")

# =============================================================================
# CONFIGURATION
# =============================================================================

println("="^80)
println("MONTE CARLO UNCERTAINTY ANALYSIS - SLEIPNER CO2 INJECTION")
println("="^80)

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
    n_realizations=50,
    start_time=0.0,
    end_time=15.0,
    num_snapshots=45,
    random_seed=1
)

println("\nConfiguration:")
println("  Realizations: $(config.n_realizations)")
println("  Time range: $(config.start_time) - $(config.end_time) years")
println("  Snapshots: $(config.num_snapshots)")
println("  Random seed: $(config.random_seed)")

# =============================================================================
# SETUP RESERVOIR GEOMETRY AND INJECTION
# =============================================================================

println("\n" * "-"^80)
println("Setting up reservoir geometry...")
println("-"^80)

# Load Sleipner topography and create domain
boundary_condition = :closed
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

println("  Domain: $(domain.nx) × $(domain.ny) × $(domain.nz)")
println("  Number of layers: $(length(layers))")

# Load injection location from feeder chimney data
xy, (utm_x, utm_y, depth) = load_feeder_location(topography)
println("  Injection location: grid indices $xy")
println("  UTM coordinates: ($utm_x, $utm_y)")
println("  Depth: $depth m")

# Generate injection schedule
injection_events = generate_sleipner_injection_events(layers, xy)
println("  Number of injection events: $(length(injection_events))")

# =============================================================================
# RUN MONTE CARLO SIMULATION
# =============================================================================

println("\n" * "-"^80)
println("Running Monte Carlo simulation...")
println("-"^80)

ensemble = run_monte_carlo_simulation(
    config,
    layers,
    domain,
    injection_events;
    verbose=true
);

println("\nMonte Carlo simulation completed!")

# =============================================================================
# UNCERTAINTY ANALYSIS
# =============================================================================

println("\n" * "-"^80)
println("Generating uncertainty analysis plots...")
println("-"^80)

# Create output directory for plots
mkpath("monte_carlo_results")

# 1. Total CO2 storage uncertainty
println("\n1. Plotting total CO2 storage uncertainty...")
plot_total_storage_uncertainty(
    ensemble;
    output_file="monte_carlo_results/total_storage_uncertainty.png",
    show_individual=false,
    figsize=(1000, 600)
)

# 2. Total CO2 leakage uncertainty
println("\n2. Plotting total CO2 leakage uncertainty...")
plot_total_leakage_uncertainty(
    ensemble;
    output_file="monte_carlo_results/total_leakage_uncertainty.png",
    show_individual=false,
    figsize=(1000, 600)
)

# 3. Layer-by-layer uncertainty
println("\n3. Plotting layer-by-layer uncertainty...")
n_layers = length(layers)
plot_all_layers_uncertainty(
    ensemble,
    n_layers;
    output_file="monte_carlo_results/all_layers_uncertainty.png",
    figsize=(1400, 1000)
)

# # 4. Individual layer plots for detailed analysis
# println("\n4. Generating individual layer plots...")
# for layer_idx in 1:n_layers
#     plot_layer_storage_uncertainty(
#         ensemble,
#         layer_idx;
#         layer_name="Layer $layer_idx",
#         output_file="monte_carlo_results/layer_$(layer_idx)_uncertainty.png",
#         show_individual=false,
#         figsize=(800, 600)
#     )
# end

# # 5. Plot with individual realizations (for first few layers as examples)
# println("\n5. Generating detailed plots with individual realizations...")
# for layer_idx in 1:min(3, n_layers)
#     plot_layer_storage_uncertainty(
#         ensemble,
#         layer_idx;
#         layer_name="Layer $layer_idx",
#         output_file="monte_carlo_results/layer_$(layer_idx)_with_realizations.png",
#         show_individual=true,
#         figsize=(800, 600)
#     )
# end

# 6. Layer distribution barplot at final timepoint
println("\n6. Generating layer distribution barplot...")
plot_layer_distribution_barplot(
    ensemble;
    output_file="monte_carlo_results/layer_distribution_barplot.png",
    figsize=(1000, 600)
)

# =============================================================================
# STATISTICAL SUMMARY
# =============================================================================

# Print comprehensive statistical summary
print_uncertainty_summary(ensemble)

# # =============================================================================
# # SAVE RESULTS
# # =============================================================================

# println("\n" * "-"^80)
# println("Saving ensemble data...")
# println("-"^80)

# # Save raw data for post-processing
# using JLD2

# @save "monte_carlo_results/ensemble_data.jld2" ensemble

# println("  Saved ensemble data to: monte_carlo_results/ensemble_data.jld2")

# # Export summary statistics to CSV for external analysis
# using CSV, DataFrames

# # Export total storage statistics
# storage_stats = compute_uncertainty_statistics(ensemble, :total_stored_m3)
# storage_df = DataFrame(
#     time=storage_stats.timepoints,
#     mean=storage_stats.mean,
#     median=storage_stats.median,
#     std=storage_stats.std,
#     p5=storage_stats.p5,
#     p25=storage_stats.p25,
#     p75=storage_stats.p75,
#     p95=storage_stats.p95,
#     min=storage_stats.min,
#     max=storage_stats.max
# )
# CSV.write("monte_carlo_results/total_storage_statistics.csv", storage_df)
# println("  Saved total storage statistics to: monte_carlo_results/total_storage_statistics.csv")

# # Export total leakage statistics
# leakage_stats = compute_uncertainty_statistics(ensemble, :total_leaked_m3)
# leakage_df = DataFrame(
#     time=leakage_stats.timepoints,
#     mean=leakage_stats.mean,
#     median=leakage_stats.median,
#     std=leakage_stats.std,
#     p5=leakage_stats.p5,
#     p25=leakage_stats.p25,
#     p75=leakage_stats.p75,
#     p95=leakage_stats.p95,
#     min=leakage_stats.min,
#     max=leakage_stats.max
# )
# CSV.write("monte_carlo_results/total_leakage_statistics.csv", leakage_df)
# println("  Saved total leakage statistics to: monte_carlo_results/total_leakage_statistics.csv")

# # Export layer-by-layer statistics
# for layer_idx in 1:n_layers
#     layer_stats = compute_layer_uncertainty_statistics(ensemble, layer_idx)
#     layer_df = DataFrame(
#         time=layer_stats.timepoints,
#         mean=layer_stats.mean,
#         median=layer_stats.median,
#         std=layer_stats.std,
#         p5=layer_stats.p5,
#         p25=layer_stats.p25,
#         p75=layer_stats.p75,
#         p95=layer_stats.p95,
#         min=layer_stats.min,
#         max=layer_stats.max
#     )
#     CSV.write("monte_carlo_results/layer_$(layer_idx)_statistics.csv", layer_df)
# end
# println("  Saved layer-by-layer statistics ($(n_layers) files)")
