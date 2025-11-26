using CairoMakie
export animate_trap_filling_multilayer, animate_trap_filling_birdseye_multilayer

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
