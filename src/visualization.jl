# Visualization functions for CO2 injection simulations

export create_layer_animation, create_ensemble_probability_animation, create_volume_timeseries_plot, create_injection_metrics_plot

"""
    create_layer_animation(topographies, texs_for_times, times, output_path;
                          title="CO2 Migration Through Layers",
                          n_cols=3, framerate=5, 
                          co2_colormap=:hot,
                          figure_size=(1600, 1300),
                          show_terrain=true)

Create an animated GIF showing CO2 migration through multiple layers.

Requires CairoMakie to be loaded in the calling scope.

# Arguments
- `topographies::Array{<:Real, 3}`: Topography data for all layers (n_layers × ny × nx)
- `texs_for_times::Vector{Vector}`: Interpolated textures from `interpolate_simulation_results`
- `times::Vector{Float64}`: Time points for animation frames
- `output_path::String`: Path where animation GIF will be saved
- `title::String`: Title for the animation
- `n_cols::Int=3`: Number of columns in grid layout
- `framerate::Int=5`: Frames per second for animation
- `co2_colormap::Symbol=:hot`: Colormap for CO2 overlay
- `figure_size::Tuple{Int,Int}=(1600, 1300)`: Figure size in pixels
- `show_terrain::Bool=true`: Whether to show terrain as contour lines

# Returns
- `fig`: The Makie Figure object (also saves animation to `output_path`)

# Example
```julia
using CairoMakie, CO2InjectionModeling

# ... run simulation ...
texs, n_filled = interpolate_simulation_results(tstructs, seqs, times_at_leakage, times)
fig = create_layer_animation(topographies, texs, times, "animation.gif")
```
"""
function create_layer_animation(
    topographies::Array{<:Real,3},
    texs_for_times::Vector,
    times::Vector{Float64},
    output_path::String;
    title::String="CO2 Migration Through Layers",
    n_cols::Int=3,
    framerate::Int=5,
    co2_colormap::Symbol=:hot,
    figure_size::Tuple{Int,Int}=(1600, 1300),
    show_terrain::Bool=true,
    inner_domain::Union{Nothing,Tuple{UnitRange{Int},UnitRange{Int}}}=nothing
)
    # Check if Makie is available
    if !isdefined(Main, :Figure) && !isdefined(Main, :CairoMakie)
        error("CairoMakie must be loaded before calling create_layer_animation. Add 'using CairoMakie' to your script.")
    end

    # Get Figure and Observable from the calling scope
    Figure = Main.Figure
    Observable = Main.Observable
    Label = Main.Label
    Axis = Main.Axis
    heatmap! = Main.heatmap!
    contour! = Main.contour!
    Colorbar = Main.Colorbar
    DataAspect = Main.DataAspect
    hidedecorations! = Main.hidedecorations!
    record = Main.record

    n_layers = size(topographies, 1)
    n_rows = ceil(Int, n_layers / n_cols)

    # Extract visualization domain
    if !isnothing(inner_domain)
        i_range, j_range = inner_domain
        vis_topographies = topographies[:, i_range, j_range]
        # Extract inner domain from texs_for_times
        vis_texs_for_times = [
            [tex[i_range, j_range] for tex in layer_texs]
            for layer_texs in texs_for_times
        ]
    else
        vis_topographies = topographies
        vis_texs_for_times = texs_for_times
    end

    # Initialize observables for each layer
    current_tex = [Observable(reverse(transpose(vis_texs_for_times[i][1]), dims=1)) for i in 1:n_layers]

    # Create figure
    fig = Figure(size=figure_size, figure_padding=20)

    # Title and time label
    Label(fig[1, 1:n_cols+1], title,
        fontsize=26, halign=:center, font=:bold, tellwidth=false)
    time_label = Label(fig[2, 1:n_cols+1], "Time: $(round(times[1], digits=2)) years",
        fontsize=20, halign=:center, tellwidth=false)

    # Create grid of subplots
    axes = []
    co2_masks = []

    # Get elevation range across all layers for consistent colorbar
    all_elevations = [vis_topographies[i, :, :] for i in 1:n_layers]
    elev_min = minimum(minimum, all_elevations)
    elev_max = maximum(maximum, all_elevations)

    for i in 1:n_layers
        row = div(i - 1, n_cols) + 3  # Start at row 3 (after title and time)
        col = mod(i - 1, n_cols) + 1

        ax = Axis(fig[row, col];
            aspect=DataAspect(),
            title="Layer $i",
            titlesize=18,
            titlegap=8)

        # Background: topography contours (if enabled)
        if show_terrain
            topo = vis_topographies[i, :, :]
            topo_reversed = reverse(transpose(topo), dims=1)
            contour!(ax, topo_reversed;
                levels=10,
                color=:black,
                linewidth=0.5,
                alpha=0.3)
        end

        # Overlay: CO2 distribution with contrasting color
        co2_mask = Observable(fill(NaN, size(current_tex[i][])))
        for idx in CartesianIndices(current_tex[i][])
            val = current_tex[i][][idx]
            co2_mask[][idx] = val > 1 ? val : NaN
        end

        heatmap!(ax, co2_mask;
            colormap=co2_colormap,
            colorrange=(1, 12),
            alpha=1.0,
            nan_color=:transparent)

        hidedecorations!(ax)
        push!(axes, ax)
        push!(co2_masks, co2_mask)
    end

    # Animation loop
    maxlength = maximum(length.(vis_texs_for_times))

    record(fig, output_path, 1:maxlength; framerate=framerate) do t
        for i in 1:n_layers
            idx = min(t, length(vis_texs_for_times[i]))
            new_tex = reverse(transpose(vis_texs_for_times[i][idx]), dims=1)
            current_tex[i][] = new_tex

            # Update CO2 mask
            for cidx in CartesianIndices(new_tex)
                val = new_tex[cidx]
                co2_masks[i][][cidx] = val > 1 ? val : NaN
            end
            Main.notify(co2_masks[i])
        end

        # Update time label
        if t <= length(times)
            time_label.text[] = "Time: $(round(times[t], digits=2)) years"
        end
    end

    println("Animation saved to $output_path")

    return fig
end

"""
    create_ensemble_probability_animation(topographies, fraction_filled, times, output_path;
                                         title="Ensemble CO2 Probability",
                                         n_cols=3, framerate=5,
                                         co2_colormap=:viridis,
                                         figure_size=(1600, 1300),
                                         show_terrain=true)

Create an animated GIF showing ensemble probability of CO2 presence through layers.
The CO2 probability is shown using a colormap from 0 (not present) to 1 (always present).

Requires CairoMakie to be loaded in the calling scope.

# Arguments
- `topographies::Array{<:Real, 3}`: Topography data for all layers (n_layers × ny × nx)
- `fraction_filled::Vector{Vector{Matrix{Float64}}}`: Probability (0-1) from `compute_ensemble_statistics`
- `times::Vector{Float64}`: Time points for animation frames
- `output_path::String`: Path where animation GIF will be saved
- `title::String`: Title for the animation
- `n_cols::Int=3`: Number of columns in grid layout
- `framerate::Int=5`: Frames per second for animation
- `co2_colormap::Symbol=:viridis`: Colormap for CO2 probability (from 0 to 1). Good options: `:viridis`, `:plasma`, `:inferno`, `:magma`
- `figure_size::Tuple{Int,Int}=(1600, 1300)`: Figure size in pixels
- `show_terrain::Bool=true`: Whether to show the underlying terrain as contour lines

# Returns
- `fig`: The Makie Figure object (also saves animation to `output_path`)

# Example
```julia
using CairoMakie, CO2InjectionModeling

# ... run ensemble simulations ...
fraction_filled, n_filled = compute_ensemble_statistics(all_texs)

# With terrain contours and viridis colormap (default)
fig = create_ensemble_probability_animation(topographies, fraction_filled, times, "ensemble.gif")

# With plasma colormap
fig = create_ensemble_probability_animation(topographies, fraction_filled, times, "ensemble_plasma.gif";
    co2_colormap=:plasma)

# Without terrain (only CO2 probability with inferno colormap)
fig = create_ensemble_probability_animation(topographies, fraction_filled, times, "ensemble_clean.gif";
    show_terrain=false, co2_colormap=:inferno)
```
"""
function create_ensemble_probability_animation(
    topographies::Array{<:Real,3},
    fraction_filled::Vector,
    times::Vector{Float64},
    output_path::String;
    title::String="Ensemble CO2 Probability",
    n_cols::Int=3,
    framerate::Int=5,
    co2_colormap::Symbol=:viridis,
    figure_size::Tuple{Int,Int}=(1600, 1300),
    show_terrain::Bool=true,
    inner_domain::Union{Nothing,Tuple{UnitRange{Int},UnitRange{Int}}}=nothing
)
    # Check if Makie is available
    if !isdefined(Main, :Figure) && !isdefined(Main, :CairoMakie)
        error("CairoMakie must be loaded before calling create_ensemble_probability_animation. Add 'using CairoMakie' to your script.")
    end

    # Get symbols from calling scope
    Figure = Main.Figure
    Observable = Main.Observable
    Label = Main.Label
    Axis = Main.Axis
    heatmap! = Main.heatmap!
    contour! = Main.contour!
    Colorbar = Main.Colorbar
    DataAspect = Main.DataAspect
    hidedecorations! = Main.hidedecorations!
    record = Main.record

    n_layers = size(topographies, 1)
    n_rows = ceil(Int, n_layers / n_cols)

    # Extract visualization domain
    if !isnothing(inner_domain)
        i_range, j_range = inner_domain
        vis_topographies = topographies[:, i_range, j_range]
        # Extract inner domain from fraction_filled
        vis_fraction_filled = [
            [frac[i_range, j_range] for frac in layer_fracs]
            for layer_fracs in fraction_filled
        ]
    else
        vis_topographies = topographies
        vis_fraction_filled = fraction_filled
    end

    # Initialize observables for each layer
    current_fraction = [Observable(reverse(transpose(vis_fraction_filled[i][1]), dims=1)) for i in 1:n_layers]

    # Create figure
    fig = Figure(size=figure_size, figure_padding=20)

    # Title and time label
    Label(fig[1, 1:n_cols+1], title,
        fontsize=26, halign=:center, font=:bold, tellwidth=false)
    time_label = Label(fig[2, 1:n_cols+1], "Time: $(round(times[1], digits=2)) years",
        fontsize=20, halign=:center, tellwidth=false)

    # Create grid of subplots
    axes = []
    co2_overlays = []
    co2_hmaps = []  # Store CO2 heatmaps for colorbar

    # Get elevation range across all layers for consistent contours
    all_elevations = [vis_topographies[i, :, :] for i in 1:n_layers]
    elev_min = minimum(minimum, all_elevations)
    elev_max = maximum(maximum, all_elevations)

    for i in 1:n_layers
        row = div(i - 1, n_cols) + 3  # Start at row 3 (after title and time)
        col = mod(i - 1, n_cols) + 1

        ax = Axis(fig[row, col];
            aspect=DataAspect(),
            title="Layer $i",
            titlesize=18,
            titlegap=8)

        # Background: topography contours (if enabled)
        if show_terrain
            topo = vis_topographies[i, :, :]
            topo_reversed = reverse(transpose(topo), dims=1)
            contour!(ax, topo_reversed;
                levels=10,
                color=:black,
                linewidth=0.5,
                alpha=0.3)
        end

        # Overlay: CO2 probability as a colormap
        # Create a matrix where 0 values are NaN (transparent) and >0 values use the colormap
        ny, nx = size(current_fraction[i][])
        co2_data = Observable(fill(NaN, ny, nx))

        # Initialize CO2 data matrix: NaN for zero, actual fraction for non-zero
        for idx in CartesianIndices(current_fraction[i][])
            frac = current_fraction[i][][idx]
            co2_data[][idx] = frac > 0 ? frac : NaN
        end

        # Plot CO2 probability with colormap (NaN values are transparent)
        co2_hm = heatmap!(ax, co2_data;
            colormap=co2_colormap,
            colorrange=(0, 1),
            nan_color=:transparent)

        hidedecorations!(ax)
        push!(axes, ax)
        push!(co2_overlays, co2_data)
        push!(co2_hmaps, co2_hm)
    end

    # Colorbar for CO2 probability
    Colorbar(fig[3:n_rows+2, n_cols+1],
        co2_hmaps[1],
        label="CO2 Probability",
        labelsize=16,
        ticklabelsize=14,
        width=20)

    # Animation loop
    maxlength = maximum(length.(vis_fraction_filled))

    record(fig, output_path, 1:maxlength; framerate=framerate) do t
        for i in 1:n_layers
            idx = min(t, length(vis_fraction_filled[i]))
            new_fraction = reverse(transpose(vis_fraction_filled[i][idx]), dims=1)
            current_fraction[i][] = new_fraction

            # Update CO2 data matrix with new fraction values (NaN for zero, fraction for non-zero)
            for cidx in CartesianIndices(new_fraction)
                frac = new_fraction[cidx]
                co2_overlays[i][][cidx] = frac > 0 ? frac : NaN
            end
            Main.notify(co2_overlays[i])
        end

        # Update time label
        if t <= length(times)
            time_label.text[] = "Time: $(round(times[t], digits=2)) years"
        end
    end

    println("Ensemble probability animation saved to $output_path")

    return fig
end

"""
    create_volume_timeseries_plot(tstructs, all_texs, common_times, cell_volume, output_path;
                                  title="CO2 Volume Over Time",
                                  figure_size=(1400, 1000),
                                  show_total=true,
                                  show_cumulative=false)

Create a plot showing CO2 volume evolution over time for each layer.

# Arguments
- `tstructs::Vector{TrapStructure}`: Trap structures for all layers
- `all_texs::Vector{Vector{Vector}}`: Textures for each simulation, layer, and time
- `common_times::Vector{Float64}`: Time points for all simulations
- `cell_volume::Float64`: Volume of a single grid cell (m³)
- `output_path::String`: Path where plot will be saved
- `title::String="CO2 Volume Over Time"`: Title for the plot
- `figure_size::Tuple{Int,Int}=(1400, 1000)`: Figure size in pixels
- `show_total::Bool=true`: Whether to show total volume across all layers
- `show_cumulative::Bool=false`: Whether to show cumulative volume plot

# Returns
- `fig`: The Makie Figure object (also saves plot to `output_path`)

# Example
```julia
using CairoMakie, CO2InjectionModeling

# Calculate cell volume
cell_volume = (length_x / ny) * (length_y / nx)

# Create volume timeseries plot
fig = create_volume_timeseries_plot(
    tstructs, all_texs, common_times, cell_volume,
    "volume_timeseries.png"
)
```
"""
function create_volume_timeseries_plot(
    tstructs::Vector,
    all_texs::Vector,
    common_times::Vector{Float64},
    cell_volume::Float64,
    output_path::String;
    title::String="CO2 Volume Over Time",
    figure_size::Tuple{Int,Int}=(1400, 1000),
    show_total::Bool=true,
    show_cumulative::Bool=false
)
    # Check if Makie is available
    if !isdefined(Main, :Figure) && !isdefined(Main, :CairoMakie)
        error("CairoMakie must be loaded before calling create_volume_timeseries_plot.")
    end

    # Get symbols from calling scope
    Figure = Main.Figure
    Axis = Main.Axis
    lines! = Main.lines!
    band! = Main.band!
    Label = Main.Label
    axislegend = Main.axislegend
    save = Main.save

    println("\nComputing CO2 volumes...")
    
    # Compute volumes for all simulations
    all_volumes = compute_co2_volumes_over_time(tstructs, all_texs, common_times, cell_volume)
    
    # Compute statistics
    mean_volumes, std_volumes, max_volumes = compute_ensemble_volume_statistics(all_volumes)
    
    n_layers = length(tstructs)
    n_simulations = length(all_texs)
    
    # Create figure
    n_plots = show_cumulative ? 2 : 1
    n_plots += show_total ? 1 : 0
    
    fig = Figure(size=figure_size, figure_padding=20)
    
    # Main title
    Label(fig[1, 1:2], title, fontsize=24, halign=:center, font=:bold, tellwidth=false)
    
    # Color palette for layers
    colors = Main.ColorSchemes.tab10.colors
    
    plot_row = 2
    
    # Plot 1: Individual layer volumes
    ax1 = Axis(fig[plot_row, 1:2];
        xlabel="Time (years)",
        ylabel="CO2 Volume (m³)",
        title="CO2 Volume by Layer (Mean ± Std)")
    
    for layer_idx in 1:n_layers
        color = colors[mod1(layer_idx, length(colors))]
        
        # Plot mean
        lines!(ax1, common_times, mean_volumes[layer_idx];
            label="Layer $layer_idx",
            color=color,
            linewidth=2)
        
        # Plot uncertainty band (mean ± std)
        lower = max.(mean_volumes[layer_idx] .- std_volumes[layer_idx], 0.0)
        upper = mean_volumes[layer_idx] .+ std_volumes[layer_idx]
        
        band!(ax1, common_times, lower, upper;
            color=(color, 0.2))
    end
    
    axislegend(ax1, position=:lt)
    plot_row += 1
    
    # Plot 2: Total volume (if requested)
    if show_total
        ax2 = Axis(fig[plot_row, 1:2];
            xlabel="Time (years)",
            ylabel="Total CO2 Volume (m³)",
            title="Total CO2 Volume Across All Layers")
        
        # Compute total volumes
        total_mean = sum(mean_volumes)
        total_std_squared = sum([s.^2 for s in std_volumes])  # Variance adds
        total_std = sqrt.(total_std_squared)
        
        lines!(ax2, common_times, total_mean;
            label="Mean total volume",
            color=:black,
            linewidth=3)
        
        lower_total = max.(total_mean .- total_std, 0.0)
        upper_total = total_mean .+ total_std
        
        band!(ax2, common_times, lower_total, upper_total;
            color=(:gray, 0.3),
            label="± 1 std")
        
        # Also plot individual simulation totals (lighter lines)
        for sim_idx in 1:min(n_simulations, 10)  # Limit to 10 for clarity
            sim_total = sum([all_volumes[sim_idx][layer_idx] for layer_idx in 1:n_layers])
            lines!(ax2, common_times, sim_total;
                color=(:gray, 0.3),
                linewidth=0.5)
        end
        
        axislegend(ax2, position=:lt)
        plot_row += 1
    end
    
    # Plot 3: Cumulative/stacked area (if requested)
    if show_cumulative
        ax3 = Axis(fig[plot_row, 1:2];
            xlabel="Time (years)",
            ylabel="Cumulative CO2 Volume (m³)",
            title="Cumulative CO2 Volume (Stacked by Layer)")
        
        # Create cumulative volumes
        cumulative_volumes = Vector{Vector{Float64}}(undef, n_layers)
        cumulative_volumes[1] = mean_volumes[1]
        
        for layer_idx in 2:n_layers
            cumulative_volumes[layer_idx] = cumulative_volumes[layer_idx-1] .+ mean_volumes[layer_idx]
        end
        
        # Plot stacked areas (from top to bottom)
        for layer_idx in n_layers:-1:1
            color = colors[mod1(layer_idx, length(colors))]
            
            lower = layer_idx > 1 ? cumulative_volumes[layer_idx-1] : zeros(length(common_times))
            upper = cumulative_volumes[layer_idx]
            
            band!(ax3, common_times, lower, upper;
                color=(color, 0.6),
                label="Layer $layer_idx")
        end
        
        axislegend(ax3, position=:lt, nbanks=2)
        plot_row += 1
    end
    
    # Add summary statistics
    total_final_mean = sum([mean_volumes[i][end] for i in 1:n_layers])
    println("\nVolume Statistics:")
    println("  Final time: $(round(common_times[end], digits=2)) years")
    println("  Total CO2 volume (mean): $(round(total_final_mean, sigdigits=4)) m³")
    for layer_idx in 1:n_layers
        println("    Layer $layer_idx: $(round(mean_volumes[layer_idx][end], sigdigits=4)) m³ " *
                "(std: $(round(std_volumes[layer_idx][end], sigdigits=3)) m³)")
    end
    
    save(output_path, fig)
    println("\nVolume timeseries plot saved to $output_path")
    
    return fig
end

"""
    create_injection_metrics_plot(common_times, all_volumes, injection_rate, co2_density, 
                                   output_path; title="CO2 Injection Metrics", figure_size=(1400, 800))

Create a plot showing injection-related metrics including storage efficiency and injected mass.

# Arguments
- `common_times::Vector{Float64}`: Time points (years)
- `all_volumes::Vector{Vector{Vector{Float64}}}`: Volumes from `compute_co2_volumes_over_time`
- `injection_rate::Float64`: Injection rate (kg/s)
- `co2_density::Float64`: CO2 density under storage conditions (kg/m³)
- `output_path::String`: Path where plot will be saved
- `title::String="CO2 Injection Metrics"`: Title for the plot
- `figure_size::Tuple{Int,Int}=(1400, 800)`: Figure size in pixels

# Returns
- `fig`: The Makie Figure object

# Example
```julia
# Typical CO2 density under storage conditions: ~700 kg/m³ (supercritical at ~1200m depth)
co2_density = 700.0  # kg/m³
injection_rate = 20000.0  # kg/s

fig = create_injection_metrics_plot(
    common_times, all_volumes, injection_rate, co2_density,
    "injection_metrics.png"
)
```
"""
function create_injection_metrics_plot(
    common_times::Vector{Float64},
    all_volumes::Vector,
    injection_rate::Float64,
    co2_density::Float64,
    output_path::String;
    title::String="CO2 Injection Metrics",
    figure_size::Tuple{Int,Int}=(1400, 800)
)
    # Check if Makie is available
    if !isdefined(Main, :Figure) && !isdefined(Main, :CairoMakie)
        error("CairoMakie must be loaded before calling create_injection_metrics_plot.")
    end

    # Get symbols from calling scope
    Figure = Main.Figure
    Axis = Main.Axis
    lines! = Main.lines!
    band! = Main.band!
    Label = Main.Label
    save = Main.save

    println("\nComputing injection metrics...")
    
    # Compute metrics
    cumulative_injected, mean_stored_fraction, std_stored_fraction = 
        compute_injection_metrics(common_times, all_volumes, injection_rate, co2_density)
    
    # Create figure
    fig = Figure(size=figure_size, figure_padding=20)
    
    # Main title
    Label(fig[1, 1:2], title, fontsize=24, halign=:center, font=:bold, tellwidth=false)
    
    # Plot 1: Cumulative injected mass vs time
    ax1 = Axis(fig[2, 1];
        xlabel="Time (years)",
        ylabel="Cumulative Injected Mass (tonnes)",
        title="CO2 Injection Rate")
    
    # Convert to tonnes for readability
    cumulative_injected_tonnes = cumulative_injected ./ 1000
    
    lines!(ax1, common_times, cumulative_injected_tonnes;
        color=:blue,
        linewidth=2,
        label="Cumulative injected")
    
    # Plot 2: Storage efficiency
    ax2 = Axis(fig[2, 2];
        xlabel="Time (years)",
        ylabel="Storage Efficiency",
        title="Fraction of Injected CO2 Stored (Mean ± Std)")
    
    lines!(ax2, common_times, mean_stored_fraction;
        color=:green,
        linewidth=2,
        label="Mean storage efficiency")
    
    lower = max.(mean_stored_fraction .- std_stored_fraction, 0.0)
    upper = min.(mean_stored_fraction .+ std_stored_fraction, 1.0)
    
    band!(ax2, common_times, lower, upper;
        color=(:green, 0.2),
        label="± 1 std")
    
    # Add reference lines for efficiency
    Main.hlines!(ax2, [1.0]; color=:gray, linestyle=:dash, linewidth=1, label="100% efficiency")
    Main.hlines!(ax2, [0.5]; color=:gray, linestyle=:dot, linewidth=1, label="50% efficiency")
    
    Main.axislegend(ax2, position=:rb)
    
    # Print summary
    println("\nInjection Metrics Summary:")
    println("  Total injected (final): $(round(cumulative_injected_tonnes[end], sigdigits=4)) tonnes")
    println("  Storage efficiency (final): $(round(mean_stored_fraction[end] * 100, digits=2))% " *
            "(± $(round(std_stored_fraction[end] * 100, digits=2))%)")
    println("  Injection rate: $injection_rate kg/s = $(round(injection_rate * 365.25 * 24 * 3600 / 1000, sigdigits=4)) tonnes/year")
    
    save(output_path, fig)
    println("\nInjection metrics plot saved to $output_path")
    
    return fig
end
