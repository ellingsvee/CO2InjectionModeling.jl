# Visualization functions for CO2 injection simulations

export create_layer_animation, create_ensemble_probability_animation

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
    show_terrain::Bool=true
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

    # Initialize observables for each layer
    current_tex = [Observable(reverse(transpose(texs_for_times[i][1]), dims=1)) for i in 1:n_layers]

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
    all_elevations = [topographies[i, :, :] for i in 1:n_layers]
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
            topo = topographies[i, :, :]
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
    maxlength = maximum(length.(texs_for_times))

    record(fig, output_path, 1:maxlength; framerate=framerate) do t
        for i in 1:n_layers
            idx = min(t, length(texs_for_times[i]))
            new_tex = reverse(transpose(texs_for_times[i][idx]), dims=1)
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
    show_terrain::Bool=true
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

    # Initialize observables for each layer
    current_fraction = [Observable(reverse(transpose(fraction_filled[i][1]), dims=1)) for i in 1:n_layers]

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
    all_elevations = [topographies[i, :, :] for i in 1:n_layers]
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
            topo = topographies[i, :, :]
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
    maxlength = maximum(length.(fraction_filled))

    record(fig, output_path, 1:maxlength; framerate=framerate) do t
        for i in 1:n_layers
            idx = min(t, length(fraction_filled[i]))
            new_fraction = reverse(transpose(fraction_filled[i][idx]), dims=1)
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
