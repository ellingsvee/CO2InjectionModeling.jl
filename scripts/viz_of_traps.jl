using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using SurfaceWaterIntegratedModeling: numtraps
using CairoMakie
using Colors
using Printf

# Setup
boundary_condition = :closed
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

# Select layer to plot
layer_idx = 1
layer = layers[layer_idx]

"""
    plot_trap_footprints(layer, domain; kwargs...)

Create a professional plot showing trap footprints for a single layer.
Traps are shown with a uniform color, with terrain contours underneath.

# Arguments
- `layer`: A Layer struct
- `domain`: Domain3D struct

# Keyword Arguments
- `output_file`: Path to save figure (default: nothing, returns figure)
- `trap_color`: Color for trap footprints (default: :royalblue)
- `trap_alpha`: Alpha transparency for trap footprints (default: 0.7)
- `show_contours`: Whether to show terrain contours (default: true)
- `contour_levels`: Number of contour lines (default: 12)
- `contour_color`: Color for contour lines (default: :gray50)
- `contour_linewidth`: Line width for contours (default: 0.8)
- `title`: Plot title (default: layer name)
- `coords_in_km`: If true, show coordinates in km (default: true)
- `transpose_layout`: If true, rotate plot 90° (default: true)
- `fontsize_title`: Font size for title (default: 14)
- `fontsize_labels`: Font size for axis labels (default: 12)
- `fontsize_ticks`: Font size for tick labels (default: 12)
- `figure_size`: Figure size as (width, height) tuple (default: (600, 400))
"""
function plot_trap_footprints(
    layer::Layer,
    domain::Domain3D;
    output_file::Union{String, Nothing} = nothing,
    trap_color = :royalblue,
    trap_alpha::Float64 = 0.7,
    show_contours::Bool = true,
    contour_levels::Int = 12,
    contour_color = :gray50,
    contour_linewidth::Float64 = 0.8,
    title::Union{String, Nothing} = nothing,
    coords_in_km::Bool = true,
    transpose_layout::Bool = true,
    fontsize_title::Int = 14,
    fontsize_labels::Int = 12,
    fontsize_ticks::Int = 12,
    figure_size::Tuple{Int, Int} = (600, 400)
)
    tstruct = layer.trap_structure
    num_traps_val = numtraps(tstruct)
    pad = layer.boundary_padding

    # Get grid size
    topo_size = size(tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Set up coordinate arrays
    if coords_in_km
        orig_x_coords = range(0, nx * domain.dx / 1000, length=nx)
        orig_y_coords = range(0, ny * domain.dy / 1000, length=ny)
        orig_x_label = "Easting (km)"
        orig_y_label = "Northing (km)"
    else
        orig_x_coords = range(0, nx * domain.dx, length=nx)
        orig_y_coords = range(0, ny * domain.dy, length=ny)
        orig_x_label = "X (m)"
        orig_y_label = "Y (m)"
    end

    # If transpose_layout, swap x and y
    if transpose_layout
        x_coords = orig_y_coords
        y_coords = orig_x_coords
        x_label = orig_y_label
        y_label = orig_x_label
    else
        x_coords = orig_x_coords
        y_coords = orig_y_coords
        x_label = orig_x_label
        y_label = orig_y_label
    end

    # Create figure
    fig = Figure(
        size = figure_size,
        backgroundcolor = :white,
        fontsize = fontsize_ticks
    )

    # plot_title = isnothing(title) ? layer.name : title
    # plot_title = nothing

    ax = Axis(fig[1, 1],
        # title = plot_title,
        titlesize = fontsize_title,
        titlefont = :bold,
        xlabel = x_label,
        ylabel = y_label,
        xlabelsize = fontsize_labels,
        ylabelsize = fontsize_labels,
        xticklabelsize = fontsize_ticks,
        yticklabelsize = fontsize_ticks,
        aspect = DataAspect(),
        xgridvisible = false,
        ygridvisible = false,
        spinewidth = 1
    )

    Makie.xlims!(ax, x_coords[1], x_coords[end])
    Makie.ylims!(ax, y_coords[1], y_coords[end])

    # Extract topography (remove padding)
    topo_padded = tstruct.topography
    topo_raw = topo_padded[pad+1:end-pad, pad+1:end-pad]

    # Transpose if needed
    if transpose_layout
        topo = permutedims(topo_raw)
    else
        topo = topo_raw
    end

    # Plot terrain contours
    if show_contours
        contour!(ax, x_coords, y_coords, topo,
            levels = contour_levels,
            color = (contour_color, 0.6),
            linewidth = contour_linewidth
        )
    end

    # Create trap footprint map
    trap_map_padded = zeros(Int, nx_padded, ny_padded)
    for trap_id in 1:num_traps_val
        footprint = tstruct.footprints[trap_id]
        for idx in footprint
            trap_map_padded[idx] = trap_id
        end
    end

    # Remove padding
    trap_map_raw = trap_map_padded[pad+1:end-pad, pad+1:end-pad]

    # Transpose if needed
    if transpose_layout
        trap_map = permutedims(trap_map_raw)
    else
        trap_map = trap_map_raw
    end

    # Create binary trap presence map (1.0 where trap exists, 0.0 elsewhere)
    trap_presence = Float64.(trap_map .> 0)

    # Plot trap footprints using heatmap with transparent/colored colormap
    heatmap!(ax, x_coords, y_coords, trap_presence,
        colormap = [RGBAf(0, 0, 0, 0), (trap_color, trap_alpha)],
        colorrange = (0.0, 1.0)
    )

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end

# Generate the trap footprint visualization for layer 1
println("Generating trap footprint visualization for $(layer.name)...")
fig = plot_trap_footprints(
    layer,
    domain;
    output_file="trap_footprints_open.svg",
    trap_color=colorant"#386624",
    trap_alpha=1.0,
    show_contours=true,
    contour_levels=12,
    contour_color=:gray50,
    contour_linewidth=0.8,
    coords_in_km=true,
    transpose_layout=true,
    fontsize_title=16,
    fontsize_labels=16,
    fontsize_ticks=16,
    # figure_size=(350, 600)
    figure_size=(600, 350)
)