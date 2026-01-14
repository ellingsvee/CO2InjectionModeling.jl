"""
Uncertainty analysis and visualization for Monte Carlo simulation results.

Provides functions for:
- Computing uncertainty statistics
- Visualizing uncertainty bands
- Analyzing parameter sensitivities
"""

using Plots
using CairoMakie
using Statistics
using Printf

include("monte_carlo_runner.jl")

"""
    UncertaintyStatistics

Statistical summary of Monte Carlo results at each timepoint.
"""
struct UncertaintyStatistics
    timepoints::Vector{Float64}
    mean::Vector{Float64}
    median::Vector{Float64}
    std::Vector{Float64}
    p5::Vector{Float64}   # 5th percentile
    p95::Vector{Float64}  # 95th percentile
    p25::Vector{Float64}  # 25th percentile
    p75::Vector{Float64}  # 75th percentile
    min::Vector{Float64}
    max::Vector{Float64}
end

"""
    compute_uncertainty_statistics(ensemble::MonteCarloEnsemble, field::Symbol)

Compute comprehensive uncertainty statistics for a specific field.

# Arguments
- `ensemble::MonteCarloEnsemble`: Monte Carlo ensemble results
- `field::Symbol`: Field to analyze (e.g., :total_stored_m3)

# Returns
- `UncertaintyStatistics`: Statistical summary
"""
function compute_uncertainty_statistics(ensemble::MonteCarloEnsemble, field::Symbol)
    data = extract_timeseries(ensemble, field)
    n_timepoints = size(data, 1)

    stats = UncertaintyStatistics(
        ensemble.timepoints,
        vec(mean(data, dims=2)),
        vec(median(data, dims=2)),
        vec(std(data, dims=2)),
        vec([quantile(data[i, :], 0.05) for i in 1:n_timepoints]),
        vec([quantile(data[i, :], 0.95) for i in 1:n_timepoints]),
        vec([quantile(data[i, :], 0.25) for i in 1:n_timepoints]),
        vec([quantile(data[i, :], 0.75) for i in 1:n_timepoints]),
        vec(minimum(data, dims=2)),
        vec(maximum(data, dims=2))
    )

    return stats
end

"""
    compute_layer_uncertainty_statistics(ensemble::MonteCarloEnsemble, layer_idx::Int)

Compute uncertainty statistics for a specific layer.

# Arguments
- `ensemble::MonteCarloEnsemble`: Monte Carlo ensemble results
- `layer_idx::Int`: Layer index

# Returns
- `UncertaintyStatistics`: Statistical summary for the layer
"""
function compute_layer_uncertainty_statistics(ensemble::MonteCarloEnsemble, layer_idx::Int)
    data = extract_layer_timeseries(ensemble, layer_idx)
    n_timepoints = size(data, 1)

    stats = UncertaintyStatistics(
        ensemble.timepoints,
        vec(mean(data, dims=2)),
        vec(median(data, dims=2)),
        vec(std(data, dims=2)),
        vec([quantile(data[i, :], 0.05) for i in 1:n_timepoints]),
        vec([quantile(data[i, :], 0.95) for i in 1:n_timepoints]),
        vec([quantile(data[i, :], 0.25) for i in 1:n_timepoints]),
        vec([quantile(data[i, :], 0.75) for i in 1:n_timepoints]),
        vec(minimum(data, dims=2)),
        vec(maximum(data, dims=2))
    )

    return stats
end

"""
    plot_uncertainty_band(
        stats::UncertaintyStatistics;
        title::String="",
        xlabel::String="Time (years)",
        ylabel::String="",
        output_file::Union{String,Nothing}=nothing,
        show_individual::Bool=false,
        individual_data::Union{Matrix{Float64},Nothing}=nothing,
        figsize::Tuple{Int,Int}=(800, 600)
    )

Plot uncertainty bands showing mean, median, and percentile ranges.

# Arguments
- `stats::UncertaintyStatistics`: Statistical summary to plot
- `title::String=""`: Plot title
- `xlabel::String="Time (years)"`: X-axis label
- `ylabel::String=""`: Y-axis label
- `output_file::Union{String,Nothing}=nothing`: Save to file if provided
- `show_individual::Bool=false`: Show individual realizations (gray lines)
- `individual_data::Union{Matrix{Float64},Nothing}=nothing`: Individual realization data
- `figsize::Tuple{Int,Int}=(800, 600)`: Figure size in pixels

# Returns
- `Plots.Plot`: The generated plot
"""
function plot_uncertainty_band(
    stats::UncertaintyStatistics;
    title::String="",
    xlabel::String="Time (years)",
    ylabel::String="",
    output_file::Union{String,Nothing}=nothing,
    show_individual::Bool=false,
    individual_data::Union{Matrix{Float64},Nothing}=nothing,
    figsize::Tuple{Int,Int}=(800, 600)
)
    p = plot(
        size=figsize,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        legend=:topleft,
        grid=true
    )

    # Plot individual realizations if requested
    if show_individual && individual_data !== nothing
        for i in 1:size(individual_data, 2)
            plot!(p, stats.timepoints, individual_data[:, i],
                  color=:gray, alpha=0.1, label="", linewidth=0.5)
        end
    end

    # Plot 90% confidence interval (5th-95th percentile)
    plot!(p, stats.timepoints, stats.p5,
          fillrange=stats.p95, fillalpha=0.2, fillcolor=:blue,
          label="90% CI", color=:blue, linewidth=0, linealpha=0)

    # Plot interquartile range (25th-75th percentile)
    plot!(p, stats.timepoints, stats.p25,
          fillrange=stats.p75, fillalpha=0.3, fillcolor=:blue,
          label="IQR (25%-75%)", color=:blue, linewidth=0, linealpha=0)

    # Plot median
    plot!(p, stats.timepoints, stats.median,
          label="Median", color=:red, linewidth=2)

    # Plot mean
    plot!(p, stats.timepoints, stats.mean,
          label="Mean", color=:black, linewidth=2, linestyle=:dash)

    # Save if output file specified
    if output_file !== nothing
        savefig(p, output_file)
        println("  Saved: $output_file")
    end

    return p
end

"""
    plot_total_storage_uncertainty(
        ensemble::MonteCarloEnsemble;
        output_file::Union{String,Nothing}="total_storage_uncertainty.png",
        show_individual::Bool=false,
        figsize::Tuple{Int,Int}=(800, 600)
    )

Plot uncertainty in total CO2 storage over time.
"""
function plot_total_storage_uncertainty(
    ensemble::MonteCarloEnsemble;
    output_file::Union{String,Nothing}="total_storage_uncertainty.png",
    show_individual::Bool=false,
    figsize::Tuple{Int,Int}=(800, 600)
)
    stats = compute_uncertainty_statistics(ensemble, :total_stored_m3)
    data = show_individual ? extract_timeseries(ensemble, :total_stored_m3) : nothing

    return plot_uncertainty_band(
        stats;
        title="Total CO₂ Storage - Uncertainty Analysis",
        ylabel="Stored CO₂ (m³)",
        output_file=output_file,
        show_individual=show_individual,
        individual_data=data,
        figsize=figsize
    )
end

"""
    plot_total_leakage_uncertainty(
        ensemble::MonteCarloEnsemble;
        output_file::Union{String,Nothing}="total_leakage_uncertainty.png",
        show_individual::Bool=false,
        figsize::Tuple{Int,Int}=(800, 600)
    )

Plot uncertainty in total CO2 leakage over time.
"""
function plot_total_leakage_uncertainty(
    ensemble::MonteCarloEnsemble;
    output_file::Union{String,Nothing}="total_leakage_uncertainty.png",
    show_individual::Bool=false,
    figsize::Tuple{Int,Int}=(800, 600)
)
    stats = compute_uncertainty_statistics(ensemble, :total_leaked_m3)
    data = show_individual ? extract_timeseries(ensemble, :total_leaked_m3) : nothing

    return plot_uncertainty_band(
        stats;
        title="Total CO₂ Leakage - Uncertainty Analysis",
        ylabel="Leaked CO₂ (m³)",
        output_file=output_file,
        show_individual=show_individual,
        individual_data=data,
        figsize=figsize
    )
end

"""
    plot_layer_storage_uncertainty(
        ensemble::MonteCarloEnsemble,
        layer_idx::Int;
        output_file::Union{String,Nothing}=nothing,
        show_individual::Bool=false,
        figsize::Tuple{Int,Int}=(800, 600)
    )

Plot uncertainty in CO2 storage for a specific layer.
"""
function plot_layer_storage_uncertainty(
    ensemble::MonteCarloEnsemble,
    layer_idx::Int;
    layer_name::String="Layer $layer_idx",
    output_file::Union{String,Nothing}=nothing,
    show_individual::Bool=false,
    figsize::Tuple{Int,Int}=(800, 600)
)
    stats = compute_layer_uncertainty_statistics(ensemble, layer_idx)
    data = show_individual ? extract_layer_timeseries(ensemble, layer_idx) : nothing

    return plot_uncertainty_band(
        stats;
        title="CO₂ Storage in $layer_name - Uncertainty Analysis",
        ylabel="Stored CO₂ (m³)",
        output_file=output_file,
        show_individual=show_individual,
        individual_data=data,
        figsize=figsize
    )
end

"""
    plot_all_layers_uncertainty(
        ensemble::MonteCarloEnsemble,
        n_layers::Int;
        output_file::Union{String,Nothing}="all_layers_uncertainty.png",
        figsize::Tuple{Int,Int}=(1200, 800)
    )

Create a multi-panel plot showing uncertainty for all layers.
"""
function plot_all_layers_uncertainty(
    ensemble::MonteCarloEnsemble,
    n_layers::Int;
    output_file::Union{String,Nothing}="all_layers_uncertainty.png",
    figsize::Tuple{Int,Int}=(1200, 800)
)
    # Create subplots
    plots = []

    for layer_idx in 1:n_layers
        stats = compute_layer_uncertainty_statistics(ensemble, layer_idx)

        p = plot(
            xlabel=(layer_idx == n_layers || layer_idx == n_layers - 1) ? "Time (years)" : "",
            ylabel="CO₂ (m³)",
            title="Layer $layer_idx",
            legend=false,
            grid=true
        )

        # Plot uncertainty bands
        plot!(p, stats.timepoints, stats.p5,
              fillrange=stats.p95, fillalpha=0.2, fillcolor=:blue,
              color=:blue, linewidth=0, linealpha=0)

        plot!(p, stats.timepoints, stats.p25,
              fillrange=stats.p75, fillalpha=0.3, fillcolor=:blue,
              color=:blue, linewidth=0, linealpha=0)

        plot!(p, stats.timepoints, stats.median,
              color=:red, linewidth=2)

        push!(plots, p)
    end

    # Combine into grid
    n_cols = 3
    n_rows = Int(ceil(n_layers / n_cols))

    combined_plot = plot(plots..., layout=(n_rows, n_cols), size=figsize)

    if output_file !== nothing
        savefig(combined_plot, output_file)
        println("  Saved: $output_file")
    end

    return combined_plot
end

"""
    print_uncertainty_summary(ensemble::MonteCarloEnsemble)

Print a comprehensive summary of uncertainty analysis results.
"""
function print_uncertainty_summary(ensemble::MonteCarloEnsemble)
    println("\n" * "="^80)
    println("MONTE CARLO UNCERTAINTY ANALYSIS SUMMARY")
    println("="^80)

    println("\nSimulation Configuration:")
    println("  Number of realizations: $(ensemble.config.n_realizations)")
    println("  Time range: $(ensemble.config.start_time) - $(ensemble.config.end_time) years")
    println("  Number of snapshots: $(ensemble.config.num_snapshots)")

    println("\nParameter Distributions:")
    for (param, dist) in ensemble.config.param_distributions
        println("  $param: $dist")
    end

    # Final time statistics
    final_time_idx = length(ensemble.timepoints)
    final_time = ensemble.timepoints[end]

    println("\n" * "-"^80)
    println("Results at Final Time (t = $final_time years)")
    println("-"^80)

    # Total storage
    storage_data = extract_timeseries(ensemble, :total_stored_m3)[:, :]
    final_storage = storage_data[final_time_idx, :]

    println("\nTotal CO₂ Storage:")
    @printf("  Mean:   %12.2e m³\n", mean(final_storage))
    @printf("  Median: %12.2e m³\n", median(final_storage))
    @printf("  Std:    %12.2e m³\n", std(final_storage))
    @printf("  Min:    %12.2e m³\n", minimum(final_storage))
    @printf("  Max:    %12.2e m³\n", maximum(final_storage))
    @printf("  P5:     %12.2e m³\n", quantile(final_storage, 0.05))
    @printf("  P95:    %12.2e m³\n", quantile(final_storage, 0.95))

    # Total leakage
    leakage_data = extract_timeseries(ensemble, :total_leaked_m3)[:, :]
    final_leakage = leakage_data[final_time_idx, :]

    println("\nTotal CO₂ Leakage:")
    @printf("  Mean:   %12.2e m³\n", mean(final_leakage))
    @printf("  Median: %12.2e m³\n", median(final_leakage))
    @printf("  Std:    %12.2e m³\n", std(final_leakage))
    @printf("  Min:    %12.2e m³\n", minimum(final_leakage))
    @printf("  Max:    %12.2e m³\n", maximum(final_leakage))
    @printf("  P5:     %12.2e m³\n", quantile(final_leakage, 0.05))
    @printf("  P95:    %12.2e m³\n", quantile(final_leakage, 0.95))

    # Layer-by-layer statistics
    n_layers = length(ensemble.results[1].snapshots[1].stored_by_layer_m3)

    println("\n" * "-"^80)
    println("CO₂ Distribution by Layer (t = $final_time years)")
    println("-"^80)
    println("\nLayer | Mean (m³)      | Median (m³)    | Std (m³)       | P5-P95 Range (m³)")
    println("-"^80)

    for layer_idx in 1:n_layers
        layer_data = extract_layer_timeseries(ensemble, layer_idx)
        final_layer = layer_data[final_time_idx, :]

        @printf("%5d | %14.2e | %14.2e | %14.2e | %14.2e - %14.2e\n",
                layer_idx,
                mean(final_layer),
                median(final_layer),
                std(final_layer),
                quantile(final_layer, 0.05),
                quantile(final_layer, 0.95))
    end

    println("="^80)
end

function plot_layer_distribution_barplot(
    ensemble::MonteCarloEnsemble;
    output_file::Union{String,Nothing}="layer_distribution_barplot.svg",
    title::Union{String,Nothing}=nothing,
    bar_color=CairoMakie.colorant"#212C53",
    bar_alpha::Float64=0.75,
    error_bar_color=:black,
    error_bar_linewidth::Float64=2.5,
    show_values::Bool=true,
    fontsize_title::Int=18,
    fontsize_labels::Int=14,
    fontsize_ticks::Int=12,
    fontsize_values::Int=10,
    figure_size::Tuple{Int,Int}=(700, 500),
    bar_width::Float64=0.7
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

        # Get stored amount in each layer
        for layer_idx in 1:n_layers
            layer_amount = ensemble.results[real_idx].snapshots[final_time_idx].stored_by_layer_m3[layer_idx]
            percentages[layer_idx, real_idx] = (layer_amount / total) * 100.0
        end
    end

    # Compute statistics for each layer
    mean_pct = vec(mean(percentages, dims=2))
    std_pct = vec(std(percentages, dims=2))
    p5_pct = vec([quantile(percentages[i, :], 0.05) for i in 1:n_layers])
    p95_pct = vec([quantile(percentages[i, :], 0.95) for i in 1:n_layers])

    # Create figure
    fig = CairoMakie.Figure(
        size=figure_size,
        backgroundcolor=:white,
        fontsize=fontsize_ticks
    )

    # Create axis
    title_row = isnothing(title) ? 1 : 2
    ax = CairoMakie.Axis(fig[title_row, 1],
        xlabel="Layer",
        ylabel="Total stored (%)",
        xlabelsize=fontsize_labels,
        ylabelsize=fontsize_labels,
        xticklabelsize=fontsize_ticks,
        yticklabelsize=fontsize_ticks,
        xticks=(1:n_layers, ["L$i" for i in 1:n_layers]),
        xgridvisible=false,
        ygridvisible=true,
        ygridcolor=(:gray, 0.3),
        ygridstyle=:dot,
        spinewidth=1.5
    )

    # Calculate error bar lengths
    lower_errors = mean_pct .- p5_pct
    upper_errors = p95_pct .- mean_pct

    # Plot bars with alpha transparency
    CairoMakie.barplot!(ax, 1:n_layers, mean_pct,
        color=(bar_color, bar_alpha),
        width=bar_width,
        strokewidth=0
    )

    # Plot error bars (90% confidence interval)
    CairoMakie.errorbars!(ax, 1:n_layers, mean_pct, lower_errors, upper_errors,
        color=error_bar_color,
        linewidth=error_bar_linewidth,
        whiskerwidth=10
    )

    # Add percentage value labels on top of bars
    if show_values
        for i in 1:n_layers
            label_y = mean_pct[i] + upper_errors[i] + maximum(mean_pct) * 0.02
            CairoMakie.text!(ax, i, label_y,
                text=@sprintf("%.1f%%", mean_pct[i]),
                align=(:center, :bottom),
                fontsize=fontsize_values,
                color=:black
            )
        end
    end

    # Add title if provided
    if !isnothing(title)
        CairoMakie.Label(fig[1, 1], title, fontsize=fontsize_title, font=:bold)
    end

    # Adjust y-axis limits to make room for labels
    if show_values
        CairoMakie.ylims!(ax, 0, maximum(mean_pct .+ upper_errors) * 1.15)
    end

    # Save if output file specified
    if !isnothing(output_file)
        CairoMakie.save(output_file, fig)
        println("  Saved: $output_file")
    end

    return fig
end

function print_percentages(
    ensemble::MonteCarloEnsemble;
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

    # Compute statistics for each layer
    mean_pct = vec(mean(percentages, dims=2))
    std_pct = vec(std(percentages, dims=2))

    println("\nLayer | Mean (%)      | Std (%)       |")
    println("-"^60)
    for layer_idx in 1:n_layers
        @printf("%5d | %12.2f | %12.2f\n",
                layer_idx,
                mean_pct[layer_idx],
                std_pct[layer_idx])
    end
end