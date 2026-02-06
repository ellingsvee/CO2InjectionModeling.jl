using Statistics
using Printf
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints, SpillEvent

export animate_multi_layer_filling, animate_multi_layer_saturation
export plot_layer_volumes_timeseries, plot_layer_fractions_timeseries
export animate_multi_layer_filling_ensemble
export plot_layer_topographies
export plot_final_co2_distribution, create_injection_locations_dict
export plot_scenario_storage_timeseries, plot_scenario_layer_distribution, ScenarioData
export plot_scenario_storage_comparison, plot_scenario_leakage_comparison


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
    animate_multi_layer_saturation(layers, seqs, leakage_states, domain; kwargs...)

Create an animated bird's eye view of CO2 saturation in all layers simultaneously.
Displays layers in a grid layout (3x3 for 9 layers).

This visualization shows CO2 saturation (volume/capacity) instead of height,
which makes residual drainage visible: after leakage starts, the saturation
in draining traps decreases from 1.0 down to the residual saturation level.

Parameters:
- `layers`: Vector of Layer structs
- `seqs`: Vector{Vector{SpillEvent}} from fill_layers
- `leakage_states`: Vector{LeakageState} from fill_layers
- `domain`: Domain3D struct
- `output_file`: Path to save animation (default: "multi_layer_saturation.gif")
- `num_frames`: Number of frames in animation (default: 30)
- `start_time`: Start time for animation (default: 0.0)
- `end_time`: End time for animation (default: auto-detect from seqs)
- `fps`: Frames per second (default: 2)
- `colormap`: Colormap for saturation (default: :viridis)
"""
function animate_multi_layer_saturation(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    leakage_states::Vector{LeakageState},
    domain::Domain3D,
    reservoir_properties::Vector{ReservoirProperties};
    output_file::String="multi_layer_saturation.gif",
    num_frames::Int=30,
    start_time::Float64=0.0,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=2,
    colormap::Symbol=:viridis
)
    n_layers = length(layers)
    @assert n_layers == length(seqs) "Number of layers must match number of sequences"
    @assert n_layers == length(leakage_states) "Number of layers must match number of leakage states"

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

        # Compute trap capacities for saturation calculation
        trap_capacities = [tstruct.trapvolumes[i] - tstruct.subvolumes[i] for i in 1:num_traps]

        push!(layer_data, (
            tstruct=tstruct,
            num_traps=num_traps,
            tstates=tstates,
            trap_capacities=trap_capacities,
            pad=layer.boundary_padding,
            name=layer.name,
            leakage_state=leakage_states[layer_idx]
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
    saturation_observables = []
    axes = []

    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        saturation_data = Observable(zeros(Float64, nx, ny))
        push!(saturation_observables, saturation_data)

        layer_name = isnothing(layer_data[layer_idx]) ? "Layer $layer_idx" : layer_data[layer_idx].name

        ax = Axis(fig[row, col],
            title=layer_name,
            xlabel=col == 1 ? "X (m)" : "",
            ylabel=row == n_rows ? "Y (m)" : "",
            aspect=DataAspect(),
            xticklabelsvisible=(row == n_rows),
            yticklabelsvisible=(col == 1))

        hm = heatmap!(ax, x_coords, y_coords, saturation_data,
            colormap=colormap,
            colorrange=(0.0, 1.0))

        push!(axes, ax)
    end

    # Add a shared colorbar
    Colorbar(fig[:, n_cols+1], colormap=colormap, colorrange=(0.0, 1.0),
        label="CO2 Saturation")

    # Add overall title with time
    time_label = Observable("Time: 0.0 years")
    Label(fig[0, :], time_label, fontsize=20, font=:bold)

    println("Creating animation...")

    # Create the animation
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]

        total_draining_all_layers = 0
        max_saturation_all_layers = 0.0

        for layer_idx in 1:n_layers
            ld = layer_data[layer_idx]

            if isnothing(ld)
                # Empty layer - just zeros
                saturation_observables[layer_idx][] = zeros(Float64, nx, ny)
                continue
            end

            filled, volumes, _ = ld.tstates[frame_idx]
            leakage_state = ld.leakage_state

            # Create saturation map (on padded grid, then remove padding)
            saturation_map_padded = zeros(Float64, nx_padded, ny_padded)

            max_saturation = 0.0
            num_draining = 0

            for trap_id in 1:ld.num_traps
                volume = volumes[trap_id]

                # For draining traps, compute drainage-adjusted volume
                if leakage_state.draining[trap_id]
                    drained_vol = compute_volume_with_drainage(trap_id, tp, leakage_state)
                    if !isnothing(drained_vol)
                        volume = drained_vol
                    end
                    num_draining += 1
                end

                # Compute saturation = volume / capacity
                capacity = ld.trap_capacities[trap_id]
                if capacity > 0.0
                    saturation = min(1.0, volume / capacity)
                else
                    saturation = 0.0
                end

                max_saturation = max(max_saturation, saturation)

                # Fill the trap footprint with this saturation
                if saturation > 0.0
                    footprint = ld.tstruct.footprints[trap_id]
                    for idx in footprint
                        saturation_map_padded[idx] = max(saturation_map_padded[idx], saturation)
                    end
                end
            end

            # Remove padding
            saturation_map = saturation_map_padded[pad+1:end-pad, pad+1:end-pad]

            # Rescale the saturation_map to account for the sand_irreducible_water_saturation
            saturation_map *= 1 - reservoir_properties[layer_idx].sand_irreducible_water_saturation

            # Update observable
            saturation_observables[layer_idx][] = saturation_map

            total_draining_all_layers += num_draining
            max_saturation_all_layers = max(max_saturation_all_layers, max_saturation)
        end

        # Update time label
        time_label[] = "Time: $(round(tp, digits=2)) years"

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
    plot_layer_fractions_timeseries(snapshots; kwargs...)

Plot the percentage of CO2 in each layer as a function of time on a single plot.
Uses a stacked area chart to show the distribution across layers.

Parameters:
- `snapshots`: Vector of ReservoirSnapshot from generate_reservoir_snapshots
- `output_file`: Path to save figure (default: nothing, returns figure)
- `title`: Plot title (default: "CO2 Distribution Across Layers")
- `show_leaked`: If true, shows leaked fraction as top area (default: true)
- `colormap`: Colormap to use for layers (default: :viridis)
"""
function plot_layer_fractions_timeseries(
    snapshots::Vector{ReservoirSnapshot};
    output_file::Union{String,Nothing}=nothing,
    title::String="CO2 Distribution Across Layers",
    show_leaked::Bool=true,
    colormap::Symbol=:viridis
)
    n_snapshots = length(snapshots)
    n_layers = length(snapshots[1].layer_snapshots)

    # Extract time series data
    times = [s.timestamp for s in snapshots]

    # Compute fractions relative to total injected (to show where all CO2 went)
    fractions_by_layer = Matrix{Float64}(undef, n_snapshots, n_layers)
    leaked_fraction = zeros(n_snapshots)

    for (t_idx, s) in enumerate(snapshots)
        total_injected = s.total_injected_m3
        if total_injected > 0
            for layer_idx in 1:n_layers
                fractions_by_layer[t_idx, layer_idx] = s.stored_by_layer_m3[layer_idx] / total_injected * 100
            end
            leaked_fraction[t_idx] = s.total_leaked_m3 / total_injected * 100
        else
            fractions_by_layer[t_idx, :] .= 0.0
            leaked_fraction[t_idx] = 0.0
        end
    end

    # Get layer names
    layer_names = [snapshots[1].layer_snapshots[i].layer_name for i in 1:n_layers]

    # Create figure
    fig = Figure(size=(800, 500))

    ax = Axis(fig[1, 1],
        title=title,
        xlabel="Time (years)",
        ylabel="Fraction of Injected CO2 (%)",
        limits=(nothing, nothing, 0, 105))

    # Generate colors from colormap
    cmap = cgrad(colormap, n_layers + (show_leaked ? 1 : 0), categorical=true)

    # Build cumulative sums for stacked area chart
    cumulative = zeros(n_snapshots)
    bands = []
    labels = String[]

    for layer_idx in 1:n_layers
        lower = copy(cumulative)
        upper = cumulative .+ fractions_by_layer[:, layer_idx]

        b = band!(ax, times, lower, upper, color=cmap[layer_idx], label=layer_names[layer_idx])
        push!(bands, b)
        push!(labels, layer_names[layer_idx])

        cumulative = upper
    end

    # Add leaked fraction on top
    if show_leaked
        lower = copy(cumulative)
        upper = cumulative .+ leaked_fraction

        b = band!(ax, times, lower, upper, color=(:red, 0.6), label="Leaked")
        push!(bands, b)
        push!(labels, "Leaked")
    end

    # Add legend
    Legend(fig[1, 2], bands, labels, "Location", framevisible=true, labelsize=11)

    # Add horizontal line at 100%
    hlines!(ax, [100], color=:black, linestyle=:dash, linewidth=1)

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
    animate_multi_layer_filling_ensemble(ensemble, layers, domain; kwargs...)

Create an animated bird's eye view of mean CO2 plume across Monte Carlo realizations.
Displays layers in a grid layout, showing the mean plume extent at each time point.

This function requires that the ensemble was created with `store_seqs=true` in
`run_monte_carlo_simulation`.

Parameters:
- `ensemble`: MonteCarloEnsemble with stored spill sequences
- `layers`: Vector of Layer structs
- `domain`: Domain3D struct
- `output_file`: Path to save animation (default: "ensemble_mean_filling.gif")
- `num_frames`: Number of frames in animation (default: 30)
- `start_time`: Start time for animation (default: 0.0)
- `end_time`: End time for animation (default: use ensemble end_time)
- `fps`: Frames per second (default: 2)
- `show_probability`: If true, color by probability (0-100%); if false, use single color (default: false)
- `plume_color`: Color for CO2 plume when show_probability=false (default: :red)
- `probability_colormap`: Colormap for probability when show_probability=true (default: :viridis)
"""
function animate_multi_layer_filling_ensemble(
    ensemble,
    layers::Vector{Layer},
    domain::Domain3D;
    output_file::String="ensemble_mean_filling.gif",
    num_frames::Int=30,
    start_time::Float64=0.0,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=2,
    show_probability::Bool=false,
    plume_color::Symbol=:red,
    probability_colormap::Symbol=:viridis
)
    n_layers = length(layers)
    n_realizations = length(ensemble.results)

    # Check that seqs are stored
    if isnothing(ensemble.results[1].seqs)
        error("Ensemble does not contain spill sequences. Run Monte Carlo simulation with store_seqs=true")
    end

    # Use ensemble end_time if not provided
    if isnothing(end_time)
        end_time = ensemble.config.end_time
    end

    # Generate timepoints for animation
    timepoints = collect(range(start_time, stop=end_time, length=num_frames))

    # Determine grid layout
    n_cols = ceil(Int, sqrt(n_layers))
    n_rows = ceil(Int, n_layers / n_cols)

    println("Computing mean trap states for $(n_layers) layers × $(num_frames) frames × $(n_realizations) realizations...")

    # Precompute data for each layer across all realizations
    layer_data = []
    for (layer_idx, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        num_traps = numtraps(tstruct)

        # Check if this layer has any data in any realization
        has_data = any(
            !isempty(ensemble.results[r].seqs[layer_idx])
            for r in 1:n_realizations
        )

        if !has_data || num_traps == 0
            push!(layer_data, nothing)
            continue
        end

        # Compute z_vol_tables for height conversion
        z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

        # Compute trap states for all realizations at all timepoints
        realization_tstates = []
        for r in 1:n_realizations
            seq = ensemble.results[r].seqs[layer_idx]
            if isempty(seq)
                push!(realization_tstates, nothing)
            else
                tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)
                push!(realization_tstates, tstates)
            end
        end

        push!(layer_data, (
            tstruct=tstruct,
            num_traps=num_traps,
            realization_tstates=realization_tstates,
            z_vol_tables=z_vol_tables,
            pad=layer.boundary_padding,
            name=layer.name
        ))
    end

    # Get grid size from first layer
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
    presence_observables = []
    axes = []

    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        presence_data = Observable(zeros(Float64, nx, ny))
        push!(presence_observables, presence_data)

        layer_name = isnothing(layer_data[layer_idx]) ? "Layer $layer_idx" : layer_data[layer_idx].name

        ax = Axis(fig[row, col],
            title=layer_name,
            xlabel=col == 1 ? "X (m)" : "",
            ylabel=row == n_rows ? "Y (m)" : "",
            aspect=DataAspect(),
            xticklabelsvisible=(row == n_rows),
            yticklabelsvisible=(col == 1))

        # Choose colormap based on show_probability
        if show_probability
            # Use colormap to show probability gradient
            hm = heatmap!(ax, x_coords, y_coords, presence_data,
                colormap=probability_colormap,
                colorrange=(0.0, 1.0))
        else
            # Use binary colormap: transparent for no CO2, plume_color for CO2 presence
            hm = heatmap!(ax, x_coords, y_coords, presence_data,
                colormap=[RGBAf(1, 1, 1, 0), plume_color],
                colorrange=(0.0, 1.0))
        end

        push!(axes, ax)
    end

    # Add colorbar if showing probability
    if show_probability
        Colorbar(fig[:, n_cols+1],
            colormap=probability_colormap,
            colorrange=(0.0, 1.0),
            label="CO2 Presence Probability")
    end

    # Add overall title with time
    mode_label = show_probability ? "Probability" : "Mean"
    time_label = Observable("Time: 0.0 years ($mode_label across $(n_realizations) realizations)")
    Label(fig[0, :], time_label, fontsize=20, font=:bold)

    println("Creating animation...")

    # Create the animation
    mode_label = show_probability ? "Probability" : "Mean"
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]

        for layer_idx in 1:n_layers
            ld = layer_data[layer_idx]

            if isnothing(ld)
                # Empty layer - just zeros
                presence_observables[layer_idx][] = zeros(Float64, nx, ny)
                continue
            end

            # Initialize presence frequency map (on padded grid)
            # This will count how many realizations have CO2 at each location
            presence_count_padded = zeros(Int, nx_padded, ny_padded)

            for r in 1:n_realizations
                tstates = ld.realization_tstates[r]

                if isnothing(tstates)
                    continue
                end

                filled, volumes, _ = tstates[frame_idx]

                # Create binary presence map for this realization
                realization_presence_map = zeros(Bool, nx_padded, ny_padded)

                for trap_id in 1:ld.num_traps
                    volume = volumes[trap_id]
                    if volume > 0.0 || filled[trap_id]
                        footprint = ld.tstruct.footprints[trap_id]
                        for idx in footprint
                            realization_presence_map[idx] = true
                        end
                    end
                end

                # Add to count
                presence_count_padded .+= realization_presence_map
            end

            # Convert count to frequency (0.0 to 1.0)
            presence_frequency_padded = presence_count_padded ./ n_realizations

            # Remove padding
            presence_map = presence_frequency_padded[pad+1:end-pad, pad+1:end-pad]

            # Update observable
            presence_observables[layer_idx][] = presence_map
        end

        # Update time label
        time_label[] = "Time: $(round(tp, digits=2)) years ($mode_label across $(n_realizations) realizations)"

        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  Frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Animation saved to: $(output_file)")
    return nothing
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

        # Get trap states at the specified timepoint
        tstates = trap_states_at_timepoints(tstruct, seq, [time]; verbose=false)

        push!(layer_data, (
            tstruct=tstruct,
            num_traps=num_traps,
            tstate=tstates[1],  # Only one timepoint
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

        # Create CO2 presence map
        filled, volumes, _ = ld.tstate
        co2_map_padded = zeros(Float64, nx_padded, ny_padded)

        for trap_id in 1:ld.num_traps
            volume = volumes[trap_id]
            if volume > 0.0 || filled[trap_id]
                footprint = ld.tstruct.footprints[trap_id]
                for idx in footprint
                    co2_map_padded[idx] = 1.0
                end
            end
        end

        # Remove padding
        co2_map_raw = co2_map_padded[pad+1:end-pad, pad+1:end-pad]

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


"""
    create_injection_locations_dict(well_locations, well_layers)

Helper function to create the injection_locations Dict from well locations and layers.

# Arguments
- `well_locations`: Vector of CartesianIndex{2} for each well
- `well_layers`: Vector of Int layer indices for each well

# Returns
- Dict{Int, Vector{CartesianIndex{2}}} mapping layer index to well locations in that layer

# Example
```julia
locations = [CartesianIndex(32, 59), CartesianIndex(45, 70)]
layers = [1, 3]
injection_locs = create_injection_locations_dict(locations, layers)
# Returns: Dict(1 => [CartesianIndex(32, 59)], 3 => [CartesianIndex(45, 70)])
```
"""
function create_injection_locations_dict(
    well_locations::Vector{<:CartesianIndex{2}},
    well_layers::Vector{Int}
)::Dict{Int,Vector{CartesianIndex{2}}}
    @assert length(well_locations) == length(well_layers) "Locations and layers must have same length"

    result = Dict{Int,Vector{CartesianIndex{2}}}()
    for (loc, layer) in zip(well_locations, well_layers)
        if !haskey(result, layer)
            result[layer] = CartesianIndex{2}[]
        end
        push!(result[layer], loc)
    end
    return result
end


# =============================================================================
# SCENARIO COMPARISON PLOTS (for optimization experiments)
# =============================================================================

"""
    ScenarioData

Data structure for storing scenario results for comparison plotting.

# Fields
- `name::String`: Display name for the scenario (e.g., "1 Well", "2 Wells")
- `timepoints::Vector{Float64}`: Time points (years)
- `storage_mt::Vector{Float64}`: Total CO2 stored at each timepoint (Mt)
- `layer_percentages::Vector{Float64}`: Percentage of CO2 in each layer at final time
"""
struct ScenarioData
    name::String
    timepoints::Vector{Float64}
    storage_mt::Vector{Float64}
    layer_percentages::Vector{Float64}
end


"""
    plot_scenario_storage_timeseries(
        scenarios::Vector{ScenarioData};
        kwargs...
    )

Create a professional time series plot comparing CO2 storage evolution across scenarios.

# Arguments
- `scenarios`: Vector of ScenarioData with storage time series

# Keyword Arguments
- `output_file`: Path to save figure (default: nothing)
- `title`: Plot title (default: nothing)
- `colors`: Vector of colors for each scenario (default: professional palette)
- `linewidth`: Line width (default: 2.5)
- `show_markers`: Show markers at data points (default: false)
- `marker_size`: Size of markers (default: 8)
- `fontsize_title`: Title font size (default: 18)
- `fontsize_labels`: Axis label font size (default: 14)
- `fontsize_ticks`: Tick label font size (default: 12)
- `fontsize_legend`: Legend font size (default: 12)
- `figure_size`: Figure size as (width, height) (default: (700, 450))
- `show_legend`: Whether to show legend (default: true)
- `legend_position`: Legend position (default: :rb for right-bottom)

# Example
```julia
scenarios = [
    ScenarioData("1 Well", times, storage1, pct1),
    ScenarioData("2 Wells", times, storage2, pct2),
    ScenarioData("3 Wells", times, storage3, pct3)
]
plot_scenario_storage_timeseries(scenarios; output_file="storage_comparison.svg")
```
"""
function plot_scenario_storage_timeseries(
    scenarios::Vector{ScenarioData};
    output_file::Union{String,Nothing}=nothing,
    title::Union{String,Nothing}=nothing,
    colors::Union{Vector,Nothing}=nothing,
    linewidth::Float64=2.5,
    show_markers::Bool=false,
    marker_size::Int=8,
    fontsize_title::Int=18,
    fontsize_labels::Int=14,
    fontsize_ticks::Int=12,
    fontsize_legend::Int=12,
    figure_size::Tuple{Int,Int}=(700, 450),
    show_legend::Bool=true,
    legend_position=:rb
)
    # Professional color palette if not provided
    if isnothing(colors)
        colors = [
            colorant"#2E86AB",  # Blue
            colorant"#A23B72",  # Magenta
            colorant"#F18F01",  # Orange
            colorant"#C73E1D",  # Red
            colorant"#3B1F2B"   # Dark
        ]
    end

    # Create figure
    fig = Figure(
        size=figure_size,
        backgroundcolor=:white,
        fontsize=fontsize_ticks
    )

    # Create axis
    title_row = isnothing(title) ? 1 : 2
    ax = Axis(fig[title_row, 1],
        xlabel="Time (years)",
        ylabel="CO₂ Stored (Mt)",
        xlabelsize=fontsize_labels,
        ylabelsize=fontsize_labels,
        xticklabelsize=fontsize_ticks,
        yticklabelsize=fontsize_ticks,
        xgridvisible=true,
        ygridvisible=true,
        xgridcolor=(:gray, 0.3),
        ygridcolor=(:gray, 0.3),
        xgridstyle=:dot,
        ygridstyle=:dot,
        spinewidth=1.5
    )

    # Plot each scenario
    for (i, scenario) in enumerate(scenarios)
        color = colors[mod1(i, length(colors))]

        lines!(ax, scenario.timepoints, scenario.storage_mt,
            color=color,
            linewidth=linewidth,
            label=scenario.name
        )

        if show_markers
            scatter!(ax, scenario.timepoints, scenario.storage_mt,
                color=color,
                markersize=marker_size
            )
        end
    end

    # Add legend
    if show_legend
        axislegend(ax, position=legend_position,
            labelsize=fontsize_legend,
            framevisible=true,
            framecolor=(:gray, 0.5),
            bgcolor=(:white, 0.9))
    end

    # Add title
    if !isnothing(title)
        Label(fig[1, 1], title, fontsize=fontsize_title, font=:bold)
    end

    # Save
    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end


"""
    plot_scenario_layer_distribution(
        scenarios::Vector{ScenarioData};
        kwargs...
    )

Create a professional grouped bar chart comparing layer distributions across scenarios.

# Arguments
- `scenarios`: Vector of ScenarioData with layer_percentages

# Keyword Arguments
- `output_file`: Path to save figure (default: nothing)
- `title`: Plot title (default: nothing)
- `colors`: Vector of colors for each scenario (default: professional palette)
- `bar_alpha`: Bar transparency (default: 0.9)
- `show_values`: Show percentage values on bars (default: true)
- `fontsize_title`: Title font size (default: 18)
- `fontsize_labels`: Axis label font size (default: 14)
- `fontsize_ticks`: Tick label font size (default: 12)
- `fontsize_values`: Value label font size (default: 9)
- `fontsize_legend`: Legend font size (default: 11)
- `figure_size`: Figure size as (width, height) (default: (800, 450))
- `bar_width`: Width of individual bars (default: 0.2)
- `legend_position`: Legend position (default: :rt)

# Example
```julia
scenarios = [
    ScenarioData("1 Well", times, storage1, [10.0, 15.0, 20.0, ...]),
    ScenarioData("2 Wells", times, storage2, [12.0, 18.0, 22.0, ...])
]
plot_scenario_layer_distribution(scenarios; output_file="layer_comparison.svg")
```
"""
function plot_scenario_layer_distribution(
    scenarios::Vector{ScenarioData};
    output_file::Union{String,Nothing}=nothing,
    title::Union{String,Nothing}=nothing,
    colors::Union{Vector,Nothing}=nothing,
    bar_alpha::Float64=0.9,
    show_values::Bool=true,
    fontsize_title::Int=18,
    fontsize_labels::Int=14,
    fontsize_ticks::Int=12,
    fontsize_values::Int=9,
    fontsize_legend::Int=11,
    figure_size::Tuple{Int,Int}=(800, 450),
    bar_width::Float64=0.2,
    legend_position=:rt
)
    # Professional color palette if not provided
    #let brown = rgb("#271B11")
    #let green = rgb("#386624")
    #let orange = rgb("#A49841")
    #let blue = rgb("#74AFB9")

    if isnothing(colors)
        colors = [
            colorant"#386624",  # Green
            colorant"#A49841",  # Orange
            colorant"#74AFB9",   # Blue
            colorant"#271B11",  # Brown
        ]
    end


    # Remove the first scenario
    scenarios_plot = scenarios[2:end]

    n_scenarios = length(scenarios_plot)
    n_layers = length(scenarios[1].layer_percentages)

    # Create figure
    fig = Figure(
        size=figure_size,
        backgroundcolor=:white,
        fontsize=fontsize_ticks
    )

    # Create axis
    title_row = isnothing(title) ? 1 : 2
    ax = Axis(fig[title_row, 1],
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


    # Plot bars for each scenario
    for (s_idx, scenario) in enumerate(scenarios_plot)
        color = colors[mod1(s_idx, length(colors))]

        # Offset for this scenario's bars within each group
        offset = (s_idx - (n_scenarios + 1) / 2) * bar_width
        positions = collect(1:n_layers) .+ offset

        barplot!(ax, positions, scenario.layer_percentages,
            color=(color, bar_alpha),
            width=bar_width * 0.9,
            strokewidth=0,
            label=scenario.name
        )

        # Add value labels
        if show_values
            for (pos, val) in zip(positions, scenario.layer_percentages)
                if val > 2.0  # Only show labels for bars with meaningful values
                    text!(ax, pos, val + 1.0,
                        text=@sprintf("%.0f", val),
                        align=(:center, :bottom),
                        fontsize=fontsize_values,
                        color=:black
                    )
                end
            end
        end
    end

    # Add legend
    axislegend(ax, position=legend_position,
        labelsize=fontsize_legend,
        framevisible=true,
        framecolor=(:gray, 0.5),
        bgcolor=(:white, 0.9))

    # Add title
    if !isnothing(title)
        Label(fig[1, 1], title, fontsize=fontsize_title, font=:bold)
    end

    # Adjust y-axis limits to accommodate labels
    if show_values
        max_val = maximum(maximum(s.layer_percentages) for s in scenarios_plot)
        ylims!(ax, 0, max_val * 1.15)
    end

    # Save
    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end


"""
    plot_scenario_storage_comparison(
        scenario_names::Vector{String},
        storage_values::Vector{Float64};
        kwargs...
    )

Create a professional bar chart comparing total CO2 storage across scenarios.

# Arguments
- `scenario_names`: Names for each scenario
- `storage_values`: Total storage in Mt for each scenario

# Keyword Arguments
- `output_file`: Path to save figure (default: nothing)
- `title`: Plot title (default: nothing)
- `colors`: Vector of colors for each bar (default: professional palette)
- `bar_alpha`: Bar transparency (default: 1.0)
- `show_values`: Show values on bars (default: true)
- `ylabel`: Y-axis label (default: "CO₂ Stored (Mt)")
- `fontsize_labels`: Axis label font size (default: 16)
- `fontsize_ticks`: Tick label font size (default: 14)
- `fontsize_values`: Value label font size (default: 14)
- `figure_size`: Figure size (default: (600, 400))
- `bar_width`: Width of bars (default: 0.6)
"""
function plot_scenario_storage_comparison(
    scenario_names::Vector{String},
    storage_values::Vector{Float64};
    output_file::Union{String,Nothing}=nothing,
    title::Union{String,Nothing}=nothing,
    colors::Union{Vector,Nothing}=nothing,
    bar_alpha::Float64=1.0,
    show_values::Bool=true,
    ylabel::String="CO₂ Stored (Mt)",
    fontsize_labels::Int=16,
    fontsize_ticks::Int=14,
    fontsize_values::Int=14,
    figure_size::Tuple{Int,Int}=(600, 400),
    bar_width::Float64=0.6
)
    if isnothing(colors)
        colors = [
            colorant"#386624",  # Green
            colorant"#A49841",  # Orange
            colorant"#74AFB9",  # Blue
            colorant"#271B11",  # Brown
        ]
    end

    n_scenarios = length(scenario_names)

    fig = Figure(
        size=figure_size,
        backgroundcolor=:white,
        fontsize=fontsize_ticks
    )

    title_row = isnothing(title) ? 1 : 2
    ax = Axis(fig[title_row, 1],
        # xlabel = "Scenario",
        ylabel=ylabel,
        xlabelsize=fontsize_labels,
        ylabelsize=fontsize_labels,
        xticklabelsize=fontsize_ticks,
        yticklabelsize=fontsize_ticks,
        xticks=(1:n_scenarios, scenario_names),
        xticklabelrotation=0.0,
        xgridvisible=false,
        ygridvisible=true,
        ygridcolor=(:gray, 0.3),
        ygridstyle=:dot,
        spinewidth=1.5
    )

    # Plot bars
    bar_colors = [(colors[mod1(i, length(colors))], bar_alpha) for i in 1:n_scenarios]
    barplot!(ax, 1:n_scenarios, storage_values,
        color=bar_colors,
        width=bar_width,
        strokewidth=0
    )

    # Add value labels
    if show_values
        for (i, val) in enumerate(storage_values)
            text!(ax, i, val + maximum(storage_values) * 0.02,
                text=@sprintf("%.2f", val),
                align=(:center, :bottom),
                fontsize=fontsize_values,
                color=:black
            )
        end
    end

    # Adjust y-axis
    if show_values
        ylims!(ax, 0, maximum(storage_values) * 1.12)
    end

    if !isnothing(title)
        Label(fig[1, 1], title, fontsize=18, font=:bold)
    end

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end


"""
    plot_scenario_leakage_comparison(
        scenario_names::Vector{String},
        leakage_values::Vector{Float64};
        kwargs...
    )

Create a professional bar chart comparing total CO2 leakage across scenarios.

# Arguments
- `scenario_names`: Names for each scenario
- `leakage_values`: Total leakage in Mt for each scenario

# Keyword Arguments
- Same as plot_scenario_storage_comparison, with ylabel default "CO₂ Leaked (Mt)"
"""
function plot_scenario_leakage_comparison(
    scenario_names::Vector{String},
    leakage_values::Vector{Float64};
    output_file::Union{String,Nothing}=nothing,
    title::Union{String,Nothing}=nothing,
    colors::Union{Vector,Nothing}=nothing,
    bar_alpha::Float64=1.0,
    show_values::Bool=true,
    ylabel::String="CO₂ Leaked (Mt)",
    fontsize_labels::Int=16,
    fontsize_ticks::Int=14,
    fontsize_values::Int=14,
    figure_size::Tuple{Int,Int}=(600, 400),
    bar_width::Float64=0.6
)
    if isnothing(colors)
        colors = [
            colorant"#386624",  # Green
            colorant"#A49841",  # Orange
            colorant"#74AFB9",  # Blue
            colorant"#271B11",  # Brown
        ]
    end

    n_scenarios = length(scenario_names)

    fig = Figure(
        size=figure_size,
        backgroundcolor=:white,
        fontsize=fontsize_ticks
    )

    title_row = isnothing(title) ? 1 : 2
    ax = Axis(fig[title_row, 1],
        # xlabel = "Scenario",
        ylabel=ylabel,
        xlabelsize=fontsize_labels,
        ylabelsize=fontsize_labels,
        xticklabelsize=fontsize_ticks,
        yticklabelsize=fontsize_ticks,
        xticks=(1:n_scenarios, scenario_names),
        xticklabelrotation=0.0,
        xgridvisible=false,
        ygridvisible=true,
        ygridcolor=(:gray, 0.3),
        ygridstyle=:dot,
        spinewidth=1.5
    )

    # Plot bars
    bar_colors = [(colors[mod1(i, length(colors))], bar_alpha) for i in 1:n_scenarios]
    barplot!(ax, 1:n_scenarios, leakage_values,
        color=bar_colors,
        width=bar_width,
        strokewidth=0
    )

    # Add value labels
    if show_values
        max_val = maximum(leakage_values)
        for (i, val) in enumerate(leakage_values)
            text!(ax, i, val + max_val * 0.02,
                text=@sprintf("%.2f", val),
                align=(:center, :bottom),
                fontsize=fontsize_values,
                color=:black
            )
        end
    end

    # Adjust y-axis
    if show_values && maximum(leakage_values) > 0
        ylims!(ax, 0, maximum(leakage_values) * 1.15)
    end

    if !isnothing(title)
        Label(fig[1, 1], title, fontsize=18, font=:bold)
    end

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end
