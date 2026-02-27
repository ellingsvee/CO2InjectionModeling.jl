module CO2BatchFillMakieExt

using CO2BatchFill
using Makie
using Statistics
using Printf
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints, SpillEvent

# Import the function stubs to extend them
import CO2BatchFill: animate_multi_layer_filling, plot_layer_volumes_timeseries
import CO2BatchFill: plot_layer_topographies, plot_final_co2_distribution

"""
    animate_multi_layer_filling(layers, seqs, domain; kwargs...)

Create an animated bird's eye view of CO2 filling in all layers simultaneously.
Displays layers in a 3x3 grid layout (for 9 layers).

Parameters:
- `layers`: Vector of Layer structs
- `seqs`: Vector{Vector{SpillEvent}} from fill_layers
- `domain`: Domain3D struct
- `output_file`: Path to save animation (default: "multi_layer_filling.gif")
- `num_frames`: Number of frames in animation (default: 30)
- `start_time`: Start time for animation (default: 0.0)
- `end_time`: End time for animation (default: auto-detect from seqs)
- `fps`: Frames per second (default: 2)
- `colormap`: Colormap for heights (default: :thermal)
- `max_CO2_height`: Maximum CO2 height for colorscale (default: 20.0)
"""
function animate_multi_layer_filling(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    domain::Domain3D;
    output_file::String="multi_layer_filling.gif",
    num_frames::Int=30,
    start_time::Float64=0.0,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=2,
    colormap::Symbol=:thermal,
    max_CO2_height::Float64=20.0
)
    n_layers = length(layers)
    @assert n_layers == length(seqs) "Number of layers must match number of sequences"

    # Determine grid layout (3x3 for 9 layers, or calculate appropriate layout)
    n_cols = ceil(Int, sqrt(n_layers))
    n_rows = ceil(Int, n_layers / n_cols)

    # Determine end_time from all sequences if not provided
    if isnothing(end_time)
        max_times = Float64[]
        for seq in seqs
            if !isempty(seq)
                max_t = maximum(se.timestamp for se in seq)
                if isfinite(max_t)
                    push!(max_times, max_t)
                elseif length(seq) > 1
                    push!(max_times, seq[end-1].timestamp + 1.0)
                end
            end
        end
        end_time = isempty(max_times) ? 15.0 : maximum(max_times)
    end

    # Generate timepoints for animation
    timepoints = collect(range(start_time, stop=end_time, length=num_frames))

    # Precompute data for each layer
    println("Computing trap states for $(n_layers) layers × $(num_frames) frames...")

    layer_data = []
    for (layer_idx, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        num_traps = numtraps(tstruct)
        seq = seqs[layer_idx]

        # Skip empty layers
        if isempty(seq) || num_traps == 0
            push!(layer_data, nothing)
            continue
        end

        # Get trap states at each timepoint
        tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)

        # Compute z_vol_tables for height conversion
        z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

        push!(layer_data, (
            tstruct=tstruct,
            num_traps=num_traps,
            tstates=tstates,
            z_vol_tables=z_vol_tables,
            pad=layer.boundary_padding,
            name=layer.name
        ))
    end

    # Get grid size from first layer (assume all layers have same size after unpadding)
    first_valid_idx = findfirst(!isnothing, layer_data)
    if isnothing(first_valid_idx)
        error("No valid layers found")
    end

    first_layer = layer_data[first_valid_idx]
    pad = first_layer.pad
    topo_size = size(first_layer.tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Set up figure with grid layout
    fig = Figure(size=(400 * n_cols, 350 * n_rows + 100))

    x_coords = range(0, nx * domain.dx, length=nx)
    y_coords = range(0, ny * domain.dy, length=ny)

    # Create observables and axes for each layer
    height_observables = []
    axes = []

    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        height_data = Observable(zeros(Float64, nx, ny))
        push!(height_observables, height_data)

        layer_name = isnothing(layer_data[layer_idx]) ? "Layer $layer_idx" : layer_data[layer_idx].name

        ax = Axis(fig[row, col],
            title=layer_name,
            xlabel=col == 1 ? "X (m)" : "",
            ylabel=row == n_rows ? "Y (m)" : "",
            aspect=DataAspect(),
            xticklabelsvisible=(row == n_rows),
            yticklabelsvisible=(col == 1))

        hm = heatmap!(ax, x_coords, y_coords, height_data,
            colormap=colormap,
            colorrange=(0.0, max_CO2_height))

        push!(axes, ax)
    end

    # Add a shared colorbar
    Colorbar(fig[:, n_cols+1], colormap=colormap, colorrange=(0.0, max_CO2_height),
        label="CO2 Height (m)")

    # Add overall title with time
    time_label = Observable("Time: 0.0 years")
    Label(fig[0, :], time_label, fontsize=20, font=:bold)

    println("Creating animation...")

    # Create the animation
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]

        total_volume_all_layers = 0.0
        max_height_all_layers = 0.0

        for layer_idx in 1:n_layers
            ld = layer_data[layer_idx]

            if isnothing(ld)
                # Empty layer - just zeros
                height_observables[layer_idx][] = zeros(Float64, nx, ny)
                continue
            end

            filled, volumes, _ = ld.tstates[frame_idx]

            # Create height map (on padded grid, then remove padding)
            height_map_padded = zeros(Float64, nx_padded, ny_padded)

            max_height = 0.0
            for trap_id in 1:ld.num_traps
                volume = volumes[trap_id]
                if volume > 0.0 || filled[trap_id]
                    z_vol_table = ld.z_vol_tables[trap_id]
                    height = volume_to_height(volume, trap_id, z_vol_table, ld.tstruct)
                    max_height = max(max_height, height)

                    footprint = ld.tstruct.footprints[trap_id]
                    for idx in footprint
                        height_map_padded[idx] = max(height_map_padded[idx], height)
                    end
                end
            end

            # Remove padding
            height_map = height_map_padded[pad+1:end-pad, pad+1:end-pad]

            # Update observable
            height_observables[layer_idx][] = height_map

            total_volume_all_layers += sum(volumes)
            max_height_all_layers = max(max_height_all_layers, max_height)
        end

        # Update time label
        time_label[] = "Time: $(round(tp, digits=2)) years | Max height: $(round(max_height_all_layers, digits=1)) m"

        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  Frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Animation saved to: $(output_file)")
    return nothing
end


"""
    plot_layer_volumes_timeseries(snapshots; kwargs...)

Plot CO2 volumes in each layer as a function of time using subplots.

Parameters:
- `snapshots`: Vector of ReservoirSnapshot from generate_reservoir_snapshots
- `output_file`: Path to save figure (default: nothing, returns figure)
- `title`: Overall title (default: "CO2 Volume by Layer")
- `colormap`: Colormap for layers (default: :tab10)
"""
function plot_layer_volumes_timeseries(
    snapshots::Vector{ReservoirSnapshot};
    output_file::Union{String,Nothing}=nothing,
    title::String="CO2 Volume by Layer",
    colormap::Symbol=:tab10
)
    n_snapshots = length(snapshots)
    n_layers = length(snapshots[1].layer_snapshots)

    # Extract time series data
    times = [s.timestamp for s in snapshots]
    volumes_by_layer = [
        [s.stored_by_layer_m3[i] / 1e6 for s in snapshots]  # Convert to M m³
        for i in 1:n_layers
    ]

    # Get layer names
    layer_names = [snapshots[1].layer_snapshots[i].layer_name for i in 1:n_layers]

    # Fixed 3x3 grid for 9 layers
    n_cols = 3
    n_rows = 3

    # Create figure with clean styling
    fig = Figure(
        size=(900, 800),
        backgroundcolor=:white
    )

    # Generate colors from colormap
    cmap = cgrad(colormap, n_layers, categorical=true)

    # Find global y-axis maximum for consistent scaling
    max_vol = maximum(maximum(v) for v in volumes_by_layer)
    y_max = max_vol * 1.15  # Add 15% padding

    # Create subplots for each layer
    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        # Determine which labels to show
        show_xlabel = row == n_rows
        show_ylabel = col == 1

        ax = Axis(fig[row, col],
            title=layer_names[layer_idx],
            titlesize=14,
            titlefont=:bold,
            xlabel=show_xlabel ? "Time (years)" : "",
            ylabel=show_ylabel ? "Volume (M m³)" : "",
            xlabelsize=12,
            ylabelsize=12,
            xticklabelsize=10,
            yticklabelsize=10,
            xticklabelsvisible=show_xlabel,
            yticklabelsvisible=show_ylabel,
            xgridvisible=true,
            ygridvisible=true,
            xgridcolor=(:gray, 0.3),
            ygridcolor=(:gray, 0.3),
            xgridstyle=:dot,
            ygridstyle=:dot,
            spinewidth=1,
            xtickwidth=1,
            ytickwidth=1)

        layer_color = cmap[layer_idx]

        # Fill area under curve first (so line is on top)
        band!(ax, times, zeros(n_snapshots), volumes_by_layer[layer_idx],
            color=(layer_color, 0.4))

        # Plot volume line
        lines!(ax, times, volumes_by_layer[layer_idx],
            color=layer_color,
            linewidth=2.5)

        # Add markers at data points
        scatter!(ax, times, volumes_by_layer[layer_idx],
            color=layer_color,
            markersize=4,
            strokewidth=0)

        ylims!(ax, 0, y_max)
        xlims!(ax, times[1], times[end])
    end

    # Add overall title
    Label(fig[0, :], title, fontsize=18, font=:bold, padding=(0, 0, 10, 0))

    # Adjust spacing
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end


"""
    plot_layer_topographies(layers, domain; kwargs...)

Create a grid of contour plots with heatmaps showing the topography (depth) of each layer.
Each subplot displays the top surface of the corresponding sand layer with its own colorbar.

Parameters:
- `layers`: Vector of Layer structs
- `domain`: Domain3D struct
- `output_file`: Path to save figure (default: nothing, returns figure)
- `title`: Overall title (default: "Layer Topographies")
- `colormap`: Colormap for depth (default: :viridis)
- `contour_levels`: Number of contour lines (default: 15)
- `show_contour_labels`: Whether to show labels on contours (default: true)
"""
function plot_layer_topographies(
    layers::Vector{Layer},
    domain::Domain3D;
    output_file::Union{String,Nothing}=nothing,
    title::String="Layer Topographies",
    colormap::Symbol=:viridis,
    contour_levels::Int=15,
    show_contour_labels::Bool=true
)
    n_layers = length(layers)

    # Determine grid layout (3x3 for 9 layers, or calculate appropriate layout)
    n_cols = ceil(Int, sqrt(n_layers))
    n_rows = ceil(Int, n_layers / n_cols)

    # Extract topography data for each layer
    layer_topos = []
    for layer in layers
        topo_padded = layer.trap_structure.topography
        pad = layer.boundary_padding

        # Remove padding
        topo = topo_padded[pad+1:end-pad, pad+1:end-pad]

        push!(layer_topos, (
            topo=topo,
            name=layer.name
        ))
    end

    # Get grid size from first layer
    nx, ny = size(layer_topos[1].topo)
    x_coords = range(0, nx * domain.dx, length=nx)
    y_coords = range(0, ny * domain.dy, length=ny)

    # Create figure with grid layout
    # Each cell needs space for axis + colorbar
    fig = Figure(size=(450 * n_cols, 350 * n_rows + 100))

    # Create subplots for each layer
    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        lt = layer_topos[layer_idx]

        # Calculate depth range for this layer only
        layer_depth_min = minimum(lt.topo)
        layer_depth_max = maximum(lt.topo)

        # Determine which labels to show
        show_xlabel = row == n_rows
        show_ylabel = col == 1

        ax = Axis(fig[row, col][1, 1],
            title=lt.name,
            titlesize=14,
            titlefont=:bold,
            xlabel=show_xlabel ? "X (m)" : "",
            ylabel=show_ylabel ? "Y (m)" : "",
            xlabelsize=12,
            ylabelsize=12,
            xticklabelsvisible=show_xlabel,
            yticklabelsvisible=show_ylabel,
            aspect=DataAspect())

        # Plot heatmap with layer-specific colorrange
        hm = heatmap!(ax, x_coords, y_coords, lt.topo,
            colormap=colormap,
            colorrange=(layer_depth_min, layer_depth_max))

        # Add contour lines
        contour!(ax, x_coords, y_coords, lt.topo,
            levels=contour_levels,
            color=:black,
            linewidth=1,
            alpha=0.6)

        # Add contour labels if requested
        if show_contour_labels
            contour!(ax, x_coords, y_coords, lt.topo,
                levels=contour_levels,
                color=:black,
                linewidth=0,
                labels=true,
                labelsize=8,
                labelfont=:regular)
        end

        # Add colorbar for this layer
        Colorbar(fig[row, col][1, 2],
            colormap=colormap,
            colorrange=(layer_depth_min, layer_depth_max),
            label="Depth (m)",
            width=15,
            ticklabelsize=10)
    end

    # Add overall title
    Label(fig[0, :], title, fontsize=18, font=:bold, padding=(0, 0, 10, 0))

    # Adjust spacing
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end

"""
    plot_final_co2_distribution(layers, seqs, domain; kwargs...)

Create a poster-ready static plot showing CO2 distribution at the final (or specified) timepoint
for all layers. Uses a uniform color for CO2 presence and optionally overlays terrain contours.

# Arguments
- `layers`: Vector of Layer structs
- `seqs`: Vector{Vector{SpillEvent}} from fill_layers
- `domain`: Domain3D struct

# Keyword Arguments
- `output_file`: Path to save figure (default: nothing, returns figure)
- `time`: Timepoint to plot (default: nothing, uses final time from sequences)
- `co2_color`: Color for CO2 presence - any Makie color specification (default: :dodgerblue)
- `co2_alpha`: Alpha transparency for CO2 (default: 0.8)
- `show_contours`: Whether to show terrain contours (default: true)
- `contour_levels`: Number of contour lines (default: 10)
- `contour_color`: Color for contour lines (default: :gray40)
- `contour_linewidth`: Line width for contours (default: 0.5)
- `contour_alpha`: Alpha for contour lines (default: 0.6)
- `title`: Overall title (default: nothing, no title)
- `show_utm_coords`: Use UTM coordinates on axes (default: true)
- `utm_origin`: Tuple (x, y) for UTM origin. Default is Sleipner coordinates (436800.0, 6468100.0).
  For other topographies, use `get_coordinate_origin(topography)` to get the appropriate origin.
- `coords_in_km`: If true, show coordinates in km relative to origin (default: false)
- `transpose_layout`: If true, rotate each subplot 90° (swap x/y axes) for wider plots (default: false)
- `fontsize_title`: Font size for main title (default: 24)
- `fontsize_layer_title`: Font size for layer titles (default: 14)
- `fontsize_labels`: Font size for axis labels (default: 12)
- `fontsize_ticks`: Font size for tick labels (default: 10)
- `figure_size`: Figure size as (width, height) tuple (default: auto-calculated)
- `row_gap`: Gap between rows in pixels (default: 5)
- `col_gap`: Gap between columns in pixels (default: 5)
- `injection_locations`: Dict mapping layer index to Vector of CartesianIndex{2} for injection well locations (default: nothing)
- `injection_marker`: Marker symbol for injection locations (default: :xcross)
- `injection_marker_color`: Color for injection markers (default: :red)
- `injection_marker_size`: Size of injection markers (default: 15)
- `injection_marker_strokewidth`: Stroke width for injection markers (default: 2.0)

# Returns
- `Figure`: The Makie figure object

# Example
```julia
# Basic usage
fig = plot_final_co2_distribution(
    layers, seqs, domain;
    output_file="final_co2.svg",
    time=15.0,
    co2_color=:royalblue,
    show_contours=true,
    transpose_layout=true
)

# With injection location markers
injection_locs = Dict(
    1 => [CartesianIndex(32, 59)],  # Well in layer 1
    3 => [CartesianIndex(45, 70)]   # Well in layer 3
)
fig = plot_final_co2_distribution(
    layers, seqs, domain;
    injection_locations=injection_locs,
    injection_marker=:xcross,
    injection_marker_color=:red
)
```
"""
function plot_final_co2_distribution(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    domain::Domain3D;
    output_file::Union{String,Nothing}=nothing,
    time::Union{Float64,Nothing}=nothing,
    co2_color=:dodgerblue,
    co2_alpha::Float64=0.8,
    show_contours::Bool=true,
    contour_levels::Int=10,
    contour_color=:gray40,
    contour_linewidth::Float64=0.5,
    contour_alpha::Float64=0.6,
    title::Union{String,Nothing}=nothing,
    show_utm_coords::Bool=true,
    utm_origin::Tuple{Float64,Float64}=(436800.0, 6468100.0),
    coords_in_km::Bool=false,
    transpose_layout::Bool=false,
    fontsize_title::Int=24,
    fontsize_layer_title::Int=14,
    fontsize_labels::Int=12,
    fontsize_ticks::Int=10,
    figure_size::Union{Tuple{Int,Int},Nothing}=nothing,
    row_gap::Int=5,
    col_gap::Int=5,
    # Injection location markers
    injection_locations::Union{Nothing,Dict{Int,Vector{CartesianIndex{2}}}}=nothing,
    injection_marker=:xcross,
    injection_marker_color=:red,
    injection_marker_size::Int=15,
    injection_marker_strokewidth::Float64=2.0
)
    n_layers = length(layers)
    @assert n_layers == length(seqs) "Number of layers must match number of sequences"

    # Determine grid layout (3x3 for 9 layers)
    n_cols = ceil(Int, sqrt(n_layers))
    n_rows = ceil(Int, n_layers / n_cols)

    # Determine timepoint to use
    if isnothing(time)
        # Find maximum finite time across all sequences
        max_times = Float64[]
        for seq in seqs
            if !isempty(seq)
                max_t = maximum(se.timestamp for se in seq)
                if isfinite(max_t)
                    push!(max_times, max_t)
                elseif length(seq) > 1
                    push!(max_times, seq[end-1].timestamp + 1.0)
                end
            end
        end
        time = isempty(max_times) ? 15.0 : maximum(max_times)
    end

    # Precompute CO2 masks for all layers at the specified timepoint
    co2_masks = create_plume_extent_masks(layers, seqs, domain, time; return_3d=false)

    # Precompute data for each layer
    layer_data = []
    for (layer_idx, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        num_traps = numtraps(tstruct)
        seq = seqs[layer_idx]

        # Skip empty layers
        if isempty(seq) || num_traps == 0
            push!(layer_data, nothing)
            continue
        end

        push!(layer_data, (
            tstruct=tstruct,
            num_traps=num_traps,
            pad=layer.boundary_padding,
            name=layer.name
        ))
    end

    # Get grid size from first valid layer
    first_valid_idx = findfirst(!isnothing, layer_data)
    if isnothing(first_valid_idx)
        error("No valid layers found")
    end

    first_layer = layer_data[first_valid_idx]
    pad = first_layer.pad
    topo_size = size(first_layer.tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Set up coordinate arrays
    if coords_in_km
        # Use local coordinates in km (relative to origin, compact display)
        orig_x_coords = range(0, nx * domain.dx / 1000, length=nx)
        orig_y_coords = range(0, ny * domain.dy / 1000, length=ny)
        orig_x_label = "Easting (km)"
        orig_y_label = "Northing (km)"
    elseif show_utm_coords
        # Use full UTM coordinates
        orig_x_coords = range(utm_origin[1], utm_origin[1] + nx * domain.dx, length=nx)
        orig_y_coords = range(utm_origin[2], utm_origin[2] + ny * domain.dy, length=ny)
        orig_x_label = "Easting (m)"
        orig_y_label = "Northing (m)"
    else
        # Use local coordinates in meters
        orig_x_coords = range(0, nx * domain.dx, length=nx)
        orig_y_coords = range(0, ny * domain.dy, length=ny)
        orig_x_label = "X (m)"
        orig_y_label = "Y (m)"
    end

    # If transpose_layout, swap x and y (rotate plots 90 degrees)
    if transpose_layout
        x_coords = orig_y_coords
        y_coords = orig_x_coords
        x_label = orig_y_label
        y_label = orig_x_label
        cell_aspect = nx / ny  # Inverted aspect ratio
    else
        x_coords = orig_x_coords
        y_coords = orig_y_coords
        x_label = orig_x_label
        y_label = orig_y_label
        cell_aspect = ny / nx
    end

    # Calculate figure size based on grid aspect ratio if not provided
    if isnothing(figure_size)
        # Base size per subplot
        subplot_width = 180
        subplot_height = round(Int, subplot_width * cell_aspect)

        title_space = isnothing(title) ? 0 : 40
        fig_width = n_cols * subplot_width + 80  # Extra space for labels
        fig_height = n_rows * subplot_height + title_space + 60

        figure_size = (fig_width, fig_height)
    end

    # Create figure with compact settings
    fig = Figure(
        size=figure_size,
        backgroundcolor=:white,
        fontsize=fontsize_ticks
    )

    # Create subplots for each layer
    for layer_idx in 1:n_layers
        # Calculate grid position (fill rows first, then columns)
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        # Adjust row for title if present
        plot_row = isnothing(title) ? row : row + 1

        ld = layer_data[layer_idx]
        layer_name = isnothing(ld) ? "Layer $layer_idx" : ld.name

        # Determine which labels to show (only on edges)
        show_xlabel = row == n_rows
        show_ylabel = col == 1

        ax = Axis(fig[plot_row, col],
            title=layer_name,
            titlesize=fontsize_layer_title,
            titlefont=:bold,
            xlabel=show_xlabel ? x_label : "",
            ylabel=show_ylabel ? y_label : "",
            xlabelsize=fontsize_labels,
            ylabelsize=fontsize_labels,
            xticklabelsize=fontsize_ticks,
            yticklabelsize=fontsize_ticks,
            xticklabelsvisible=show_xlabel,
            yticklabelsvisible=show_ylabel,
            aspect=DataAspect(),
            xgridvisible=false,
            ygridvisible=false,
            spinewidth=1
        )

        # Set axis limits explicitly
        xlims!(ax, x_coords[1], x_coords[end])
        ylims!(ax, y_coords[1], y_coords[end])

        if isnothing(ld)
            # Empty layer - show blank background
            continue
        end

        # Extract topography for contours (remove padding)
        topo_padded = ld.tstruct.topography
        topo_raw = topo_padded[pad+1:end-pad, pad+1:end-pad]

        # Get precomputed CO2 mask for this layer (already unpadded, convert to Float64)
        co2_map_raw = Float64.(co2_masks[layer_idx])

        # Transpose data if needed (rotate 90 degrees)
        if transpose_layout
            topo = permutedims(topo_raw)
            co2_map = permutedims(co2_map_raw)
        else
            topo = topo_raw
            co2_map = co2_map_raw
        end

        # Plot terrain contours if requested
        if show_contours
            contour!(ax, x_coords, y_coords, topo,
                levels=contour_levels,
                color=(contour_color, contour_alpha),
                linewidth=contour_linewidth
            )
        end

        # Plot CO2 presence as a heatmap with custom colormap
        # Use transparent for no CO2, and the specified color with alpha for CO2 presence
        heatmap!(ax, x_coords, y_coords, co2_map,
            colormap=[RGBAf(0, 0, 0, 0), (co2_color, co2_alpha)],
            colorrange=(0.0, 1.0)
        )

        # Plot injection location markers if provided for this layer
        if !isnothing(injection_locations) && haskey(injection_locations, layer_idx)
            locs = injection_locations[layer_idx]
            for loc in locs
                # Convert grid index to coordinates
                # Grid indices are 1-based, need to map to coordinate system
                if coords_in_km
                    loc_x = (loc[1] - 1) * domain.dx / 1000
                    loc_y = (loc[2] - 1) * domain.dy / 1000
                elseif show_utm_coords
                    loc_x = utm_origin[1] + (loc[1] - 1) * domain.dx
                    loc_y = utm_origin[2] + (loc[2] - 1) * domain.dy
                else
                    loc_x = (loc[1] - 1) * domain.dx
                    loc_y = (loc[2] - 1) * domain.dy
                end

                # Swap coordinates if transposed
                if transpose_layout
                    loc_x, loc_y = loc_y, loc_x
                end

                scatter!(ax, [loc_x], [loc_y],
                    marker=injection_marker,
                    markersize=injection_marker_size,
                    color=injection_marker_color,
                    strokewidth=injection_marker_strokewidth,
                    strokecolor=injection_marker_color
                )
            end
        end
    end

    # Add overall title if provided
    if !isnothing(title)
        Label(fig[1, :], title, fontsize=fontsize_title, font=:bold)
    end

    # Adjust spacing for compact appearance
    colgap!(fig.layout, col_gap)
    rowgap!(fig.layout, row_gap)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end

end # module
