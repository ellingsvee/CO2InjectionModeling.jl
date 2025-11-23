using CairoMakie
export plot_cross_section, animate_trap_filling, animate_trap_filling_birdseye
export animate_trap_filling_multilayer, animate_trap_filling_birdseye_multilayer

"""
Plot a 2D cross-section with lithology and/or values (e.g., CO2 saturation, pressure).

Parameters:
- lithology: 3D array (nx × ny × nz) with rock types (optional if mode=:values)
- values: 3D array (nx × ny × nz) with values to plot (optional if mode=:lithology)
- topography: SleipnerTopography struct
- dz: vertical grid spacing (meters)
- direction: :x or :y (which direction to slice)
- slice_index: index along the slice direction (default: middle)
- mode: :lithology (only lithology), :values (only values), :both (overlay), or :overlay (values with lithology contours)
- value_label: label for the value colorbar (default: "Values")
- value_colormap: colormap for values (default: :viridis)
- value_range: tuple (min, max) for value colorrange (default: auto from data)
- alpha: transparency for overlay mode (default: 0.7)

Returns:
- fig: Makie Figure object

Example:
```julia
# Plot CO2 saturation over lithology
topography = load_sleipner_topography()
lithology = reconstruct_3d_lithology(topography, 1.0)
co2_saturation = zeros(size(lithology))
# ... set some values in co2_saturation using masks ...

fig = plot_cross_section(lithology, co2_saturation, topography, 1.0,
                        direction=:x, mode=:overlay,
                        value_label="CO2 Saturation", value_colormap=:thermal)
```
"""
function plot_cross_section(
    lithology::Union{Array{Int,3}, Nothing},
    values::Union{Array{<:Real,3}, Nothing},
    domain::Domain3D;
    direction::Symbol = :x,
    slice_index::Union{Int,Nothing} = nothing,
    mode::Symbol = :both,
    value_label::String = "Values",
    value_colormap::Symbol = :viridis,
    value_range::Union{Tuple{Real,Real}, Nothing} = nothing,
    alpha::Float64 = 0.7
)
    # Validate inputs
    if mode == :lithology && lithology === nothing
        error("lithology must be provided when mode=:lithology")
    end
    if (mode == :values || mode == :both || mode == :overlay) && values === nothing
        error("values must be provided when mode=:values, :both, or :overlay")
    end
    if mode == :both || mode == :overlay
        if lithology === nothing
            error("lithology must be provided when mode=:both or :overlay")
        end
    end

    # Get dimensions from whichever array is provided
    if lithology !== nothing
        nx, ny, nz = size(lithology)
    else
        nx, ny, nz = size(values)
    end

    # Determine slice index if not provided (use middle)
    if isnothing(slice_index)
        slice_index = direction == :x ? div(nx, 2) : div(ny, 2)
    end

    # Extract 2D slices
    if direction == :x
        litho_slice = lithology !== nothing ? lithology[slice_index, :, :] : nothing
        value_slice = values !== nothing ? values[slice_index, :, :] : nothing
        xlabel_text = "Y position (m)"
        title_text = "Cross-section at x = $(round((slice_index-1)*domain.dx, digits=1)) m"
        horizontal_coords = range(0, domain.ny * domain.dy,
                                 length = litho_slice !== nothing ? size(litho_slice, 1) : size(value_slice, 1))
    elseif direction == :y
        litho_slice = lithology !== nothing ? lithology[:, slice_index, :] : nothing
        value_slice = values !== nothing ? values[:, slice_index, :] : nothing
        xlabel_text = "X position (m)"
        title_text = "Cross-section at y = $(round((slice_index-1)*domain.dy, digits=1)) m"
        horizontal_coords = range(0, domain.nx * domain.dx,
                                 length = litho_slice !== nothing ? size(litho_slice, 1) : size(value_slice, 1))
    else
        error("direction must be :x or :y")
    end

    # Create depth coordinates
    depth_coords = domain.depth_max .- (0.5:nz) .* domain.dz

    # Check if value slice has any non-zero values and compute safe color range (if applicable)
    safe_value_range = nothing
    if value_slice !== nothing
        vmin, vmax = minimum(value_slice), maximum(value_slice)
        if vmin == vmax
            # All values are identical - create a small range around the value
            @warn "The selected slice has constant values ($vmin) for '$value_label'. " *
                  "The plot will be created with a small range around this value."
            safe_value_range = (vmin - 0.5, vmin + 0.5)
        else
            safe_value_range = (vmin, vmax)
        end

        # Additional warning for all-zero slices
        if vmin == 0 && vmax == 0
            @warn "The selected slice has no non-zero values for '$value_label'."
        end
    end

    # Create figure
    fig = Figure(size = (1200, 600))
    ax = Axis(fig[1, 1],
              xlabel = xlabel_text,
              ylabel = "Depth (m)",
              title = title_text,
              yreversed = true)

    # Plot based on mode
    if mode == :lithology
        # Only lithology
        litho_colormap = [:brown, :gold, :gray]  # caprock, sand, shale
        hm = heatmap!(ax, horizontal_coords, depth_coords, Float64.(litho_slice),
                     colormap = litho_colormap,
                     colorrange = (-0.5, 2.5))
        Colorbar(fig[1, 2], hm,
                label = "Lithology",
                ticks = (0:2, ["Caprock", "Sand", "Shale"]))

    elseif mode == :values
        # Only values
        vrange = value_range !== nothing ? value_range : safe_value_range
        hm = heatmap!(ax, horizontal_coords, depth_coords, Float64.(value_slice),
                     colormap = value_colormap,
                     colorrange = vrange)
        Colorbar(fig[1, 2], hm, label = value_label)

    elseif mode == :both
        # Side-by-side plots
        fig = Figure(size = (1600, 600))

        # Lithology on left
        ax1 = Axis(fig[1, 1],
                  xlabel = xlabel_text,
                  ylabel = "Depth (m)",
                  title = title_text * " - Lithology",
                  yreversed = true)
        litho_colormap = [:brown, :gold, :gray]
        hm1 = heatmap!(ax1, horizontal_coords, depth_coords, Float64.(litho_slice),
                      colormap = litho_colormap,
                      colorrange = (-0.5, 2.5))
        Colorbar(fig[1, 2], hm1,
                label = "Lithology",
                ticks = (0:2, ["Caprock", "Sand", "Shale"]))

        # Values on right
        ax2 = Axis(fig[1, 3],
                  xlabel = xlabel_text,
                  ylabel = "Depth (m)",
                  title = title_text * " - $(value_label)",
                  yreversed = true)
        vrange = value_range !== nothing ? value_range : safe_value_range
        hm2 = heatmap!(ax2, horizontal_coords, depth_coords, Float64.(value_slice),
                      colormap = value_colormap,
                      colorrange = vrange)
        Colorbar(fig[1, 4], hm2, label = value_label)

    elseif mode == :overlay
        # Values overlaid on lithology with transparency
        # First plot lithology as background
        litho_colormap = [:brown, :gold, :gray]
        heatmap!(ax, horizontal_coords, depth_coords, Float64.(litho_slice),
                colormap = litho_colormap,
                colorrange = (-0.5, 2.5),
                alpha = 0.3)

        # Then overlay values with transparency
        vrange = value_range !== nothing ? value_range : safe_value_range
        hm_values = heatmap!(ax, horizontal_coords, depth_coords, Float64.(value_slice),
                            colormap = value_colormap,
                            colorrange = vrange,
                            alpha = alpha)

        # Add colorbar for values only
        Colorbar(fig[1, 2], hm_values, label = value_label)
    else
        error("mode must be :lithology, :values, :both, or :overlay")
    end

    return fig
end

"""
Animate trap filling over time from simulation snapshots.

This function creates an animated cross-section showing how traps fill with CO2 over time.
The masks are computed on-the-fly for each frame to avoid storing large 3D arrays.

Parameters:
- snapshots: Vector of SimulationLayerSnapshot containing spill events
- layer: Layer struct containing trap_structure
- lithology: 3D array (nx × ny × nz) with rock types
- domain: Domain3D struct
- output_file: Path where to save the animation (default: "trap_filling_animation.gif")
- direction: :x or :y (which direction to slice, default: :x)
- slice_index: index along the slice direction (default: middle)
- fps: frames per second for the animation (default: 2)
- value_colormap: colormap for CO2 saturation (default: :thermal)
- alpha: transparency for overlay mode (default: 0.7)

Returns:
- Nothing (saves animation to file)

Example:
```julia
snapshots = fill_layer_co2(layers[1].trap_structure, domain, reservoir_props, injection)
animate_trap_filling(snapshots, layers[1], lithology, domain,
                    output_file="trap_animation.gif", direction=:x)
```
"""
function animate_trap_filling(
    snapshots::Vector{SimulationLayerSnapshot},
    layer::Layer,
    lithology::Array{Int,3},
    domain::Domain3D;
    output_file::String = "trap_filling_animation.gif",
    direction::Symbol = :x,
    slice_index::Union{Int,Nothing} = nothing,
    fps::Int = 2,
    value_colormap::Symbol = :thermal,
    alpha::Float64 = 0.7
)
    nx, ny, nz = size(lithology)

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

    println("Creating animation with $(length(snapshots)) frames...")

    # Create the animation
    record(fig, output_file, eachindex(snapshots); framerate=fps) do frame_idx
        snapshot = snapshots[frame_idx]

        # Use the height matrix stored in the snapshot to create accurate 3D mask
        height_matrix = snapshot.heights

        # Create 3D mask from heights (accounts for actual water levels, not just filled traps)
        mask_3d = create_co2_mask_3d_from_heights(height_matrix, layer, domain)

        # Create CO2 saturation field
        co2_saturation_3d = zeros(Float64, nx, ny, nz)
        co2_saturation_3d[mask_3d] .= 0.6

        # Extract the slice
        co2_slice = direction == :x ? co2_saturation_3d[slice_index, :, :] : co2_saturation_3d[:, slice_index, :]

        # Get max height for display
        max_height = maximum(height_matrix)

        # Update the observable
        co2_data[] = co2_slice
        time_text[] = "$(title_base) - Time: $(round(snapshot.timestamp, digits=2)) - Max height: $(round(max_height, digits=1))m"

        if frame_idx % 10 == 0 || frame_idx == length(snapshots)
            println("  Frame $(frame_idx)/$(length(snapshots))")
        end
    end

    println("Animation saved to: $(output_file)")

    return nothing
end

"""
Animate trap filling from a bird's eye view (top-down view of the layer).

This function creates an animated 2D heatmap showing the layer from above,
with traps colored as they fill with CO2 over time.

Parameters:
- snapshots: Vector of SimulationLayerSnapshot containing spill events
- layer: Layer struct containing trap_structure
- domain: Domain3D struct
- output_file: Path where to save the animation (default: "trap_filling_birdseye.gif")
- fps: frames per second for the animation (default: 2)
- value_colormap: colormap for CO2 saturation (default: :thermal)

Returns:
- Nothing (saves animation to file)

Example:
```julia
snapshots = fill_layer_co2(layers[1].trap_structure, domain, reservoir_props, injection)
animate_trap_filling_birdseye(snapshots, layers[1], domain,
                              output_file="birdseye_animation.gif")
```
"""
function animate_trap_filling_birdseye(
    snapshots::Vector{SimulationLayerSnapshot},
    layer::Layer,
    domain::Domain3D;
    output_file::String = "trap_filling_birdseye.gif",
    fps::Int = 2,
    value_colormap::Symbol = :thermal
)
    # Get the topography dimensions from the trap structure
    trap_topo = layer.trap_structure.topography
    nx_topo, ny_topo = size(trap_topo)

    # Create coordinate ranges for the heatmap
    x_coords = range(0, domain.nx * domain.dx, length=nx_topo)
    y_coords = range(0, domain.ny * domain.dy, length=ny_topo)

    # Set up the figure and axis
    fig = Figure(size = (1200, 1000))

    # Create observable for the CO2 saturation data
    co2_data = Observable(zeros(Float64, nx_topo, ny_topo))
    time_text = Observable("Bird's Eye View - Time: 0.0")

    ax = Axis(fig[1, 1],
              xlabel = "X position (m)",
              ylabel = "Y position (m)",
              title = time_text,
              aspect = DataAspect())

    # Create the heatmap
    hm = heatmap!(ax, x_coords, y_coords, co2_data,
                  colormap = value_colormap,
                  colorrange = (0.0, 0.6))

    Colorbar(fig[1, 2], hm, label = "CO2 Saturation")

    println("Creating bird's eye view animation with $(length(snapshots)) frames...")

    # Create the animation
    record(fig, output_file, eachindex(snapshots); framerate=fps) do frame_idx
        snapshot = snapshots[frame_idx]

        # Use the height matrix stored in the snapshot
        # The height matrix is already 2D and shows CO2 heights at each location
        height_matrix = snapshot.heights

        # Create a 2D map of CO2 saturation based on height
        # Non-zero heights indicate presence of CO2
        co2_map = zeros(Float64, nx_topo, ny_topo)
        co2_map[height_matrix .> 0.0] .= 0.6

        # Get max height for display
        max_height = maximum(height_matrix)
        num_filled_locations = count(height_matrix .> 0.0)

        # Update the observable
        co2_data[] = co2_map
        time_text[] = "Bird's Eye View - Time: $(round(snapshot.timestamp, digits=2)) - Max height: $(round(max_height, digits=1))m - Cells: $(num_filled_locations)"

        if frame_idx % 10 == 0 || frame_idx == length(snapshots)
            println("  Frame $(frame_idx)/$(length(snapshots))")
        end
    end

    println("Bird's eye view animation saved to: $(output_file)")

    return nothing
end

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

        trap_topo = layers[layer_idx].trap_structure.topography
        nx_topo, ny_topo = size(trap_topo)

        x_coords = range(0, domain.nx * domain.dx, length=nx_topo)
        y_coords = range(0, domain.ny * domain.dy, length=ny_topo)

        co2_observables[layer_idx] = Observable(zeros(Float64, nx_topo, ny_topo))

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

            trap_topo = layers[layer_idx].trap_structure.topography
            nx_topo, ny_topo = size(trap_topo)

            # Create a 2D map of CO2 saturation based on height
            co2_map = zeros(Float64, nx_topo, ny_topo)
            co2_map[height_matrix .> 0.0] .= 0.6

            # Update stats
            total_filled_cells += count(height_matrix .> 0.0)
            max_height_all = max(max_height_all, maximum(height_matrix))

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
