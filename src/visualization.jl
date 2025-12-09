using CairoMakie
using Statistics
export animate_trap_filling_multilayer, animate_trap_filling_birdseye_multilayer
export plot_total_co2_volume, plot_layer_co2_volumes

"""
Animate trap filling over time for a multi-layer system.

This function creates an animated cross-section showing how traps fill with CO2
across multiple layers over time. All layers are shown in a single combined view.

Parameters:
- snapshots: Vector of SimulationSnapshot containing layer snapshots
- layers: Vector of Layer structs containing trap_structures
- lithology: 3D array (nx × ny × nz) with rock types
- domain: Domain3D struct
- output_file: Path where to save the animation (default: "trap_filling_multilayer.gif")
- direction: :x or :y (which direction to slice, default: :x)
- slice_index: index along the slice direction (default: middle)
- fps: frames per second for the animation (default: 2)
- value_colormap: colormap for CO2 saturation (default: :thermal)
- alpha: transparency for overlay mode (default: 0.7)

Returns:
- Nothing (saves animation to file)
"""
function animate_trap_filling_multilayer(
    snapshots::Vector{SimulationSnapshot},
    layers::Vector{Layer},
    lithology::Array{Int,3},
    domain::Domain3D;
    output_file::String = "trap_filling_multilayer.gif",
    direction::Symbol = :x,
    slice_index::Union{Int,Nothing} = nothing,
    fps::Int = 2,
    value_colormap::Symbol = :thermal,
    alpha::Float64 = 0.7
)
    nx, ny, nz = size(lithology)
    n_layers = length(layers)

    # Determine slice index if not provided (use middle)
    if isnothing(slice_index)
        slice_index = direction == :x ? div(nx, 2) : div(ny, 2)
    end

    # Set up the figure and axis
    fig = Figure(size = (1200, 600))

    if direction == :x
        xlabel_text = "Y position (m)"
        horizontal_coords = range(0, domain.ny * domain.dy, length=ny)
        title_base = "Cross-section at x = $(round((slice_index-1)*domain.dx, digits=1)) m"
    else
        xlabel_text = "X position (m)"
        horizontal_coords = range(0, domain.nx * domain.dx, length=nx)
        title_base = "Cross-section at y = $(round((slice_index-1)*domain.dy, digits=1)) m"
    end

    depth_coords = domain.depth_max .- (0.5:nz) .* domain.dz

    # Extract lithology slice once (it doesn't change)
    litho_slice = direction == :x ? lithology[slice_index, :, :] : lithology[:, slice_index, :]

    # Create observable for the CO2 saturation data
    co2_data = Observable(zeros(Float64, size(litho_slice)))
    time_text = Observable("Time: 0.0")

    ax = Axis(fig[1, 1],
              xlabel = xlabel_text,
              ylabel = "Depth (m)",
              title = time_text,
              yreversed = true)

    # Plot lithology as background
    litho_colormap = [:brown, :gold, :gray]
    heatmap!(ax, horizontal_coords, depth_coords, Float64.(litho_slice),
            colormap = litho_colormap,
            colorrange = (-0.5, 2.5),
            alpha = 0.3)

    # Overlay CO2 saturation
    hm_values = heatmap!(ax, horizontal_coords, depth_coords, co2_data,
                        colormap = value_colormap,
                        colorrange = (0.0, 0.6),
                        alpha = alpha)

    Colorbar(fig[1, 2], hm_values, label = "CO2 Saturation")

    println("Creating multi-layer animation with $(length(snapshots)) frames...")

    # Create the animation
    record(fig, output_file, eachindex(snapshots); framerate=fps) do frame_idx
        snapshot = snapshots[frame_idx]

        # Initialize combined CO2 saturation field
        co2_saturation_3d = zeros(Float64, nx, ny, nz)

        # Process each layer
        for layer_idx in 1:n_layers
            layer_snapshot = snapshot.layer_snapshots[layer_idx]
            height_matrix = layer_snapshot.heights

            # Create 3D mask from heights for this layer
            mask_3d = create_co2_mask_3d_from_heights(height_matrix, layers[layer_idx], domain)

            # Add CO2 saturation for this layer (use max to handle overlaps)
            co2_saturation_3d[mask_3d] .= 0.6
        end

        # Extract the slice
        co2_slice = direction == :x ? co2_saturation_3d[slice_index, :, :] : co2_saturation_3d[:, slice_index, :]

        # Update the observable
        co2_data[] = co2_slice
        time_text[] = "$(title_base) - Time: $(round(snapshot.timestamp, digits=2)) - Total CO2: $(round(snapshot.total_co2_volume, digits=2))"

        if frame_idx % 10 == 0 || frame_idx == length(snapshots)
            println("  Frame $(frame_idx)/$(length(snapshots))")
        end
    end

    println("Multi-layer animation saved to: $(output_file)")

    return nothing
end

"""
Animate trap filling from a bird's eye view for a multi-layer system.

This function creates an animated grid of 2D heatmaps showing each layer from above,
with traps colored as they fill with CO2 over time. Layers are arranged in a 3x3 grid.

Parameters:
- snapshots: Vector of SimulationSnapshot containing layer snapshots
- layers: Vector of Layer structs containing trap_structures
- domain: Domain3D struct
- output_file: Path where to save the animation (default: "trap_filling_birdseye_multilayer.gif")
- fps: frames per second for the animation (default: 2)
- value_colormap: colormap for CO2 saturation (default: :thermal)

Returns:
- Nothing (saves animation to file)
"""
function animate_trap_filling_birdseye_multilayer(
    snapshots::Vector{SimulationSnapshot},
    layers::Vector{Layer},
    domain::Domain3D;
    output_file::String = "trap_filling_birdseye_multilayer.gif",
    fps::Int = 2,
    value_colormap::Symbol = :thermal
)
    n_layers = length(layers)

    # Determine grid layout (up to 3x3 = 9 layers)
    n_cols = min(3, n_layers)
    n_rows = ceil(Int, n_layers / n_cols)

    # Set up the figure with subplots
    fig = Figure(size = (400 * n_cols, 350 * n_rows + 100))

    # Create a shared title at the top
    time_text = Observable("Bird's Eye View - Time: 0.0")
    Label(fig[0, 1:n_cols], time_text, fontsize=20)

    # Create axes and observables for each layer
    axes = Vector{Axis}(undef, n_layers)
    co2_observables = Vector{Observable{Matrix{Float64}}}(undef, n_layers)

    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        # Use original domain dimensions (without padding) for visualization
        x_coords = range(0, domain.nx * domain.dx, length=domain.nx)
        y_coords = range(0, domain.ny * domain.dy, length=domain.ny)

        co2_observables[layer_idx] = Observable(zeros(Float64, domain.nx, domain.ny))

        axes[layer_idx] = Axis(fig[row, col],
                               xlabel = "X (m)",
                               ylabel = "Y (m)",
                               title = "Layer $(layer_idx): $(layers[layer_idx].name)",
                               aspect = DataAspect())

        heatmap!(axes[layer_idx], x_coords, y_coords, co2_observables[layer_idx],
                 colormap = value_colormap,
                 colorrange = (0.0, 0.6))
    end

    # Add a shared colorbar
    Colorbar(fig[1:n_rows, n_cols+1], colormap=value_colormap, colorrange=(0.0, 0.6),
             label="CO2 Saturation")

    println("Creating bird's eye view multi-layer animation with $(length(snapshots)) frames...")

    # Create the animation
    record(fig, output_file, eachindex(snapshots); framerate=fps) do frame_idx
        snapshot = snapshots[frame_idx]

        total_filled_cells = 0
        max_height_all = 0.0

        for layer_idx in 1:n_layers
            layer_snapshot = snapshot.layer_snapshots[layer_idx]
            height_matrix = layer_snapshot.heights
            pad = layers[layer_idx].boundary_padding

            # Remove padding from height matrix for visualization
            unpadded_height = remove_boundary_padding(height_matrix, pad)

            # Create a 2D map of CO2 saturation based on height
            co2_map = zeros(Float64, domain.nx, domain.ny)

            # For cells with CO2, show presence (use 0.6 for consistency with colormap)
            has_co2 = unpadded_height .> 0.0
            co2_map[has_co2] .= 0.6

            # Update stats (count on unpadded data)
            num_filled = count(has_co2)
            total_filled_cells += num_filled
            max_height_layer = maximum(unpadded_height)
            max_height_all = max(max_height_all, max_height_layer)

            # Debug: print layer info for first and last frames
            if frame_idx == 1 || frame_idx == length(snapshots)
                if num_filled > 0
                    println("  Layer $layer_idx ($(layers[layer_idx].name)): $num_filled cells filled, max height: $(round(max_height_layer, digits=2))m, CO2 volume: $(round(layer_snapshot.co2_volume, digits=2)) m³")
                end
            end

            # Update the observable
            co2_observables[layer_idx][] = co2_map
        end

        time_text[] = "Bird's Eye View - Time: $(round(snapshot.timestamp, digits=2)) - Max height: $(round(max_height_all, digits=1))m - Total cells: $(total_filled_cells)"

        if frame_idx % 10 == 0 || frame_idx == length(snapshots)
            println("  Frame $(frame_idx)/$(length(snapshots))")
        end
    end

    println("Bird's eye view multi-layer animation saved to: $(output_file)")

    return nothing
end

"""
    plot_total_co2_volume(summary::SimulationSummary; kwargs...)
    plot_total_co2_volume(summaries::Vector{SimulationSummary}; kwargs...)

Plot the total CO2 volume over time.

For a single summary, plots the time series of total CO2 volume.
For multiple summaries, plots the mean with a shaded region showing ± one standard deviation.

# Arguments
- `summary::SimulationSummary` or `summaries::Vector{SimulationSummary}`: Simulation summary data
- `output_file::Union{String,Nothing}`: Path to save the figure (optional, displays if nothing)
- `title::String`: Plot title (default: "Total CO2 Volume Over Time")
- `xlabel::String`: X-axis label (default: "Time (years)")
- `ylabel::String`: Y-axis label (default: "CO2 Volume (m³)")
- `linecolor`: Color for the line plot (default: :blue)
- `bandcolor`: Color for the std dev band (default: (:blue, 0.3))
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 500))

# Returns
- `Figure`: The Makie figure object

# Examples
```julia
# Single simulation
fig = plot_total_co2_volume(summary)

# Multiple simulations with custom styling
fig = plot_total_co2_volume(summaries, title="CO2 Storage", linecolor=:red)

# Save to file
plot_total_co2_volume(summary, output_file="co2_volume.png")
```
"""
function plot_total_co2_volume(
    summary::SimulationSummary;
    output_file::Union{String,Nothing} = nothing,
    title::String = "Total CO2 Volume Over Time",
    xlabel::String = "Time (years)",
    ylabel::String = "CO2 Volume (m³)",
    linecolor = :blue,
    bandcolor = (:blue, 0.3),
    figsize::Tuple{Int,Int} = (800, 500)
)
    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              xlabel = xlabel,
              ylabel = ylabel,
              title = title)

    lines!(ax, summary.timepoints, summary.total_co2_volumes,
           color = linecolor, linewidth = 2)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end

function plot_total_co2_volume(
    summaries::Vector{SimulationSummary};
    output_file::Union{String,Nothing} = nothing,
    title::String = "Total CO2 Volume Over Time",
    xlabel::String = "Time (years)",
    ylabel::String = "CO2 Volume (m³)",
    linecolor = :blue,
    bandcolor = (:blue, 0.3),
    figsize::Tuple{Int,Int} = (800, 500)
)
    # Ensure all summaries have the same timepoints
    timepoints = summaries[1].timepoints
    n_points = length(timepoints)
    n_sims = length(summaries)

    # Collect volumes from all simulations
    volumes_matrix = zeros(Float64, n_points, n_sims)
    for (i, summary) in enumerate(summaries)
        volumes_matrix[:, i] = summary.total_co2_volumes
    end

    # Compute mean and std dev
    mean_volumes = mean(volumes_matrix, dims=2)[:, 1]
    std_volumes = std(volumes_matrix, dims=2)[:, 1]

    # Create figure
    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              xlabel = xlabel,
              ylabel = ylabel,
              title = title)

    # Plot shaded region for ± 1 std dev
    band!(ax, timepoints,
          mean_volumes .- std_volumes,
          mean_volumes .+ std_volumes,
          color = bandcolor)

    # Plot mean line
    lines!(ax, timepoints, mean_volumes,
           color = linecolor, linewidth = 2, label = "Mean")

    # Add legend
    axislegend(ax, position = :lt)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end

"""
    plot_layer_co2_volumes(summary::SimulationSummary; kwargs...)
    plot_layer_co2_volumes(summaries::Vector{SimulationSummary}; kwargs...)

Plot CO2 volumes in different layers over time.

For a single summary, plots time series for each layer.
For multiple summaries, plots mean values with shaded regions showing ± one standard deviation.

# Arguments
- `summary::SimulationSummary` or `summaries::Vector{SimulationSummary}`: Simulation summary data
- `output_file::Union{String,Nothing}`: Path to save the figure (optional, displays if nothing)
- `title::String`: Plot title (default: "Layer CO2 Volumes Over Time")
- `xlabel::String`: X-axis label (default: "Time (years)")
- `ylabel::String`: Y-axis label (default: "CO2 Volume (m³)")
- `layer_names::Union{Vector{String},Nothing}`: Custom layer names (default: "Layer 1", "Layer 2", etc.)
- `colors`: Color scheme for layers (default: automatic cycling)
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 500))

# Returns
- `Figure`: The Makie figure object

# Examples
```julia
# Single simulation
fig = plot_layer_co2_volumes(summary)

# Multiple simulations with custom layer names
fig = plot_layer_co2_volumes(summaries, layer_names=["Utsira", "Shale", "Lower"])

# Save to file
plot_layer_co2_volumes(summary, output_file="layer_volumes.png")
```
"""
function plot_layer_co2_volumes(
    summary::SimulationSummary;
    output_file::Union{String,Nothing} = nothing,
    title::String = "Layer CO2 Volumes Over Time",
    xlabel::String = "Time (years)",
    ylabel::String = "CO2 Volume (m³)",
    layer_names::Union{Vector{String},Nothing} = nothing,
    colors = Makie.wong_colors(),
    figsize::Tuple{Int,Int} = (800, 500)
)
    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              xlabel = xlabel,
              ylabel = ylabel,
              title = title)

    # Default layer names if not provided
    if isnothing(layer_names)
        layer_names = ["Layer $i" for i in 1:summary.num_layers]
    end

    # Plot each layer
    for layer_idx in 1:summary.num_layers
        color = colors[mod1(layer_idx, length(colors))]
        lines!(ax, summary.timepoints, summary.layer_co2_volumes[:, layer_idx],
               color = color, linewidth = 2, label = layer_names[layer_idx])
    end

    # Add legend
    axislegend(ax, position = :lt)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end

function plot_layer_co2_volumes(
    summaries::Vector{SimulationSummary};
    output_file::Union{String,Nothing} = nothing,
    title::String = "Layer CO2 Volumes Over Time",
    xlabel::String = "Time (years)",
    ylabel::String = "CO2 Volume (m³)",
    layer_names::Union{Vector{String},Nothing} = nothing,
    colors = Makie.wong_colors(),
    figsize::Tuple{Int,Int} = (800, 500)
)
    # Ensure all summaries have the same structure
    timepoints = summaries[1].timepoints
    n_points = length(timepoints)
    n_sims = length(summaries)
    n_layers = summaries[1].num_layers

    # Default layer names if not provided
    if isnothing(layer_names)
        layer_names = ["Layer $i" for i in 1:n_layers]
    end

    # Create figure
    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              xlabel = xlabel,
              ylabel = ylabel,
              title = title)

    # Process each layer
    for layer_idx in 1:n_layers
        # Collect volumes from all simulations for this layer
        layer_volumes_matrix = zeros(Float64, n_points, n_sims)
        for (i, summary) in enumerate(summaries)
            layer_volumes_matrix[:, i] = summary.layer_co2_volumes[:, layer_idx]
        end

        # Compute mean and std dev
        mean_volumes = mean(layer_volumes_matrix, dims=2)[:, 1]
        std_volumes = std(layer_volumes_matrix, dims=2)[:, 1]

        # Select color
        color = colors[mod1(layer_idx, length(colors))]
        band_color = (color, 0.3)

        # Plot shaded region for ± 1 std dev
        band!(ax, timepoints,
              mean_volumes .- std_volumes,
              mean_volumes .+ std_volumes,
              color = band_color)

        # Plot mean line
        lines!(ax, timepoints, mean_volumes,
               color = color, linewidth = 2, label = layer_names[layer_idx])
    end

    # Add legend
    axislegend(ax, position = :lt)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end
