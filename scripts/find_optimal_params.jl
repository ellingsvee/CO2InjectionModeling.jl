using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using Printf
using Plots

include("monte_carlo_analysis.jl")


# Optimizider
layer_percentages_true = [9, 10, 7, 6, 25, 8, 7, 11, 17] ./ 100

# =============================================================================
# CONFIGURATION
# =============================================================================


 config = create_default_monte_carlo_config(
    n_realizations=1000,
    start_time=0.0,
    end_time=15.0,
    num_snapshots=30,
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



function get_ensemble_with_best_fit_to_percentages(
    ensemble::MonteCarloEnsemble;
    layer_percentages_true::Vector{Float64}
)
    # Get number of layers
    n_layers = length(ensemble.results[1].snapshots[1].stored_by_layer_m3)
    n_realizations = ensemble.config.n_realizations
    final_time_idx = length(ensemble.timepoints)

    # Extract data for all layers at final timepoint
    # percentages[layer, realization]
    percentages = zeros(n_layers, n_realizations)

    for real_idx in 1:n_realizations
        # Get total stored CO2 for this realization at final time
        total = ensemble.results[real_idx].snapshots[final_time_idx].total_stored_m3
        # total = ensemble.results[real_idx].snapshots[final_time_idx].total_injected_m3
        

        # Get stored amount in each layer
        for layer_idx in 1:n_layers
            layer_amount = ensemble.results[real_idx].snapshots[final_time_idx].stored_by_layer_m3[layer_idx]
            percentages[layer_idx, real_idx] = (layer_amount / total) * 100.0
        end
    end

    # Find the realization that best fits the true percentages
    best_real_idx = 0
    best_error = Inf
    for real_idx in 1:n_realizations
        error = sum((percentages[:, real_idx] .- layer_percentages_true .* 100.0).^ 2)
        # error = sum(abs.(percentages[:, real_idx] .- layer_percentages_true .* 100.0))
        if error < best_error
            best_error = error
            best_real_idx = real_idx
        end
    end
    println("Best fitting realization: $best_real_idx with error $best_error")
    return ensemble.results[best_real_idx], best_real_idx
end

best_fit_realization, best_fit_idx = get_ensemble_with_best_fit_to_percentages(
    ensemble;
    layer_percentages_true=layer_percentages_true
);


for i in 1:9
    println(best_fit_realization.parameters[Symbol("shale_pressure_threshold_layer_$i")])
end

function plot_layer_distribution_barplot_best_fit_realization(
    ensemble::MonteCarloEnsemble,
    best_fit_idx;
    output_file::Union{String,Nothing}="layer_distribution_barplot_best_fit_realization.png",
    figsize::Tuple{Int,Int}=(1000, 600)
)
    # Get number of layers
    n_layers = length(ensemble.results[1].snapshots[1].stored_by_layer_m3)
    n_realizations = ensemble.config.n_realizations
    final_time_idx = length(ensemble.timepoints)

    # Extract data for all layers at final timepoint
    # percentages[layer, realization]
    percentages = zeros(n_layers, 1)

    # for real_idx in 1:n_realizations
    #     # Get total stored CO2 for this realization at final time
    total = ensemble.results[best_fit_idx].snapshots[final_time_idx].total_stored_m3
    # total = ensemble.results[real_idx].snapshots[final_time_idx].total_injected_m3
    

    # Get stored amount in each layer
    for layer_idx in 1:n_layers
        layer_amount = ensemble.results[best_fit_idx].snapshots[final_time_idx].stored_by_layer_m3[layer_idx]
        percentages[layer_idx, 1] = (layer_amount / total) * 100.0
    end
    # end

    # Compute statistics for each layer
    mean_pct = vec(mean(percentages, dims=2))
    std_pct = vec(std(percentages, dims=2))
    p5_pct = vec([quantile(percentages[i, :], 0.05) for i in 1:n_layers])
    p95_pct = vec([quantile(percentages[i, :], 0.95) for i in 1:n_layers])

    # Create barplot
    layer_labels = ["L$i" for i in 1:n_layers]

    p = bar(
        1:n_layers,
        mean_pct,
        # yerror=(mean_pct .- p5_pct, p95_pct .- mean_pct),
        xlabel="Layer",
        ylabel="Percentage of Total CO₂ (%)",
        title="CO₂ Distribution Across Layers at t = $(ensemble.timepoints[end]) years",
        label="Mean ± 90% CI",
        legend=:topright,
        color=:steelblue,
        alpha=0.7,
        xticks=(1:n_layers, layer_labels),
        size=figsize,
        grid=true,
        bar_width=0.6
    )

    # Add percentage labels on top of bars
    for i in 1:n_layers
        annotate!(p, i, mean_pct[i] + (p95_pct[i] - mean_pct[i]) + 1.5,
                 text(@sprintf("%.1f%%", mean_pct[i]), :center, 8))
    end

    # Save if output file specified
    if output_file !== nothing
        savefig(p, output_file)
        println("  Saved: $output_file")
    end

    return p
end

plot_layer_distribution_barplot_best_fit_realization(
    ensemble,
    best_fit_idx;
    output_file="monte_carlo_results/layer_distribution_barplot_best_fit_realization.png",
    figsize=(1000, 600)
)