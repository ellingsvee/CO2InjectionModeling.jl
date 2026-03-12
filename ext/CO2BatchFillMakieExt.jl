module CO2BatchFillMakieExt

using CO2BatchFill
using Makie
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints,
    SpillEvent, _compute_z_vol_tables
using Statistics: mean, std

import CO2BatchFill: animate_layer_filling, animate_multi_layer_filling,
    plot_layer, plot_multi_layer, plot_multi_layer_ensemble,
    plot_layer_volumes_timeseries, plot_multi_layer_volumes_timeseries,
    plot_multi_layer_ensemble_timeseries

# Shared style constants
const _CMAP = :thermal
const _MAXH = 20.0    # default max CO2 height (m)
const _LW = 4       # default linewidth
const _COLORS = (
    stored=:black,
    drained=:skyblue,
)
const _CBAR_SIZE = 25

# Private helpers
"""
Build an unpadded 2-D CO2 height array from per-trap volumes.

Optional `leaking` / `leakage_heights` keyword arguments enable physically-correct
height reporting for residually-trapped CO2 (see `height_map` in utils.jl).
"""
function _height_map(tstruct, z_vol_tables, volumes, pad, nx_p, ny_p;
    leaking::Union{Nothing,AbstractVector{Bool}}=nothing,
    leakage_heights::Union{Nothing,AbstractVector{Float64}}=nothing)
    num_traps = numtraps(tstruct)
    height_map_padded = zeros(Float64, nx_p, ny_p)
    use_leakage_correction = !isnothing(leaking) && !isnothing(leakage_heights)
    for trap_id in 1:num_traps
        vol = volumes[trap_id]
        vol <= 0.0 && continue
        h = if use_leakage_correction && leaking[trap_id] && isfinite(leakage_heights[trap_id])
            leakage_heights[trap_id]
        else
            volume_to_height(vol, trap_id, z_vol_tables[trap_id], tstruct)
        end
        h <= 0.0 && continue
        min_topo = get_min_topography_elevation(trap_id, tstruct)
        water_level = min_topo + h
        for idx in tstruct.footprints[trap_id]
            cell_h = max(0.0, water_level - tstruct.topography[idx])
            height_map_padded[idx] = max(height_map_padded[idx], cell_h)
        end
    end
    return height_map_padded[pad+1:end-pad, pad+1:end-pad]
end

"""
Draw topography contours (major + minor) onto `ax`.
"""
function _topo_contours!(ax, x_coords, y_coords, topo;
    contour_levels=10, major_contour_every=5, contour_opacity=0.8, show_labels=false)
    if contour_levels isa Integer
        mn, mx = extrema(topo)
        all_levels = collect(range(mn, mx; length=contour_levels))
    else
        all_levels = collect(contour_levels)
    end

    major_levels = all_levels[2:major_contour_every:end]
    minor_levels = setdiff(all_levels, major_levels)

    contour!(ax, x_coords, y_coords, topo;
        levels=minor_levels,
        color=(:black, contour_opacity),
        linewidth=0.8,
        labels=false
    )
    contour!(ax, x_coords, y_coords, topo;
        levels=major_levels,
        color=(:black, contour_opacity),
        linewidth=2.0,
        labels=show_labels,
        labelsize=12,
    )
end

"""
Draw a CO2 height heatmap + topography contours for one `LayerSnapshot` onto `ax`.
Returns the heatmap object.
"""
function _layer_panel!(ax, layer, snap, domain; colormap=_CMAP, max_co2_height=_MAXH,
    show_contours=true, show_labels=false, contour_levels=10, major_contour_every=5, contour_opacity=0.8,
    show_extents=false, extent_color=(:dodgerblue, 0.5))

    tstruct = layer.trap_structure
    pad = layer.pad_width
    nx_padded, ny_padded = get_padded_grid_size(layer)
    nx, ny = get_grid_size(layer)

    z_vol_tables = _compute_z_vol_tables(tstruct)
    height_array = _height_map(tstruct, z_vol_tables, snap.trap_volumes, pad, nx_padded, ny_padded;
        leaking=snap.trap_leaking,
        leakage_heights=snap.trap_leakage_height)

    x_coords = range(0.0, nx * domain.dx, length=nx)
    y_coords = range(0.0, ny * domain.dy, length=ny)
    topo = tstruct.topography[pad+1:end-pad, pad+1:end-pad]

    hm = if show_extents
        extent_array = Float64.(height_array .> 0.0)
        heatmap!(ax, x_coords, y_coords, extent_array;
            colormap=[(:white, 0.0), extent_color],
            colorrange=(0.0, 1.0),
        )
    else
        heatmap!(ax, x_coords, y_coords, height_array;
            colormap=colormap,
            colorrange=(0.0, max_co2_height),
        )
    end

    if show_contours
        _topo_contours!(ax, x_coords, y_coords, topo;
            contour_levels, major_contour_every, contour_opacity, show_labels)
    end
    return hm
end

"""
Draw stored/drained/passthrough/→next-layer volume lines for a vector of
`LayerSnapshot`s onto `ax`.
"""
function _timeseries_panel!(ax, snaps::Vector{LayerSnapshot};
    show_injected=false, linewidth=_LW, vol_scale=1.0)

    times = [s.timestamp for s in snaps]

    lines!(ax, times, [s.total_stored * vol_scale for s in snaps];
        label="Stored", color=:black, linewidth=linewidth, linestyle=:solid)
    lines!(ax, times, [s.total_drained * vol_scale for s in snaps];
        label="Drained", color=:black, linewidth=linewidth, linestyle=:dash)
end



"""
Plot CO2 height distribution for a single layer at a given snapshot time.
"""
function plot_layer(
    layer::CO2BatchFill.Layer,
    snap::CO2BatchFill.LayerSnapshot,
    domain::CO2BatchFill.Domain3D;
    output_file::Union{String,Nothing}=nothing,
    colormap::Symbol=_CMAP,
    max_co2_height::Float64=_MAXH,
    show_contours::Bool=true,
    contour_levels::Int=10,
    show_extents::Bool=false,
    extent_color=(:dodgerblue, 0.5),
)
    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1];
        xlabel="x",
        ylabel="y",
        title=snap.layer_name,
        aspect=DataAspect(),
    )
    _layer_panel!(ax, layer, snap, domain;
        colormap, max_co2_height, show_contours, contour_levels,
        show_extents, extent_color)
    if !show_extents
        Colorbar(fig[1, 2]; colormap, colorrange=(0.0, max_co2_height), label="Column height", size=_CBAR_SIZE)
    end

    if output_file !== nothing
        save(output_file, fig)
        println("Saved layer plot to: $(output_file)")
    end
    return fig
end

"""
Helper function to convert a linear index `i` into row/column indices for a given `figure_layout` (nrows, ncols).
"""
function i_to_layout(figure_layout, i)
    row = (i - 1) ÷ figure_layout[2] + 1
    col = (i - 1) % figure_layout[2] + 1
    return row, col
end

"""
Plot CO2 height distribution for multiple layers side-by-side with a shared colorbar.
"""
function plot_multi_layer(
    layers::Vector{CO2BatchFill.Layer},
    snap::CO2BatchFill.MultiLayerSnapshot,
    domain::CO2BatchFill.Domain3D;
    output_file::Union{String,Nothing}=nothing,
    colormap::Symbol=_CMAP,
    max_co2_height::Float64=_MAXH,
    show_contours::Bool=true,
    show_labels::Bool=false,
    contour_levels::Int=10,
    major_contour_every::Int=5,
    contour_opacity::Float64=0.8,
    injection_locations::Union{Nothing,Vector{Tuple{Float64,Float64}}}=nothing,
    show_leakage_locations::Bool=false,
    figure_size::Union{Tuple{Int,Int},Nothing}=nothing,
    figure_layout::Union{Tuple{Int,Int},Nothing}=nothing,
    cbar_location::Union{Tuple{Int,Int},Nothing}=nothing,
    show_extents::Bool=false,
    extent_color=(:dodgerblue, 0.5),
    colorbar_label="Column height",
)
    n = length(layers)

    figure_size = isnothing(figure_size) ? (700 * n, 600) : figure_size
    figure_layout = isnothing(figure_layout) ? (1, length(layers)) : figure_layout

    fig = Figure(size=figure_size)

    has_injection = !isnothing(injection_locations)
    has_any_leakage = false

    for i in 1:n
        figure_row, figure_col = i_to_layout(figure_layout, i)
        ax = Axis(fig[figure_row, figure_col];
            xlabel="x",
            ylabel=(i == 1 ? "y" : ""),
            title=snap.layers[i].layer_name,
            aspect=DataAspect(),
        )
        _layer_panel!(ax, layers[i], snap.layers[i], domain;
            colormap, max_co2_height, show_contours, show_labels, contour_levels, major_contour_every, contour_opacity,
            show_extents, extent_color)

        if has_injection && i == 1
            scatter!(ax, [loc[1] for loc in injection_locations], [loc[2] for loc in injection_locations];
                color=:black, marker=:xcross, markersize=35)
        end

        if show_leakage_locations
            tstruct = layers[i].trap_structure
            pad = layers[i].pad_width
            layer_snap = snap.layers[i]
            leak_xs, leak_ys = Float64[], Float64[]
            for trap_id in 1:numtraps(tstruct)
                layer_snap.trap_leaking[trap_id] || continue
                ci = find_leakage_location(trap_id, tstruct)
                x = (ci[1] - 1 - pad) * domain.dx
                y = (ci[2] - 1 - pad) * domain.dy
                push!(leak_xs, x)
                push!(leak_ys, y)
            end
            if !isempty(leak_xs)
                scatter!(ax, leak_xs, leak_ys;
                    color=:black, marker=:utriangle, markersize=30)
                has_any_leakage = true
            end
        end
    end
    if !show_extents
        cbar_location = isnothing(cbar_location) ? (1, n + 1) : cbar_location
        Colorbar(fig[cbar_location[1], cbar_location[2]]; colormap, colorrange=(0.0, max_co2_height), label=colorbar_label, size=_CBAR_SIZE)
    end

    if output_file !== nothing
        save(output_file, fig)
        println("Saved multi-layer plot to: $(output_file)")
    end
    return fig
end

# Time-series plots
"""
Plot CO2 volume time-series for a single layer.

- `vol_scale`: multiply raw SWIM volumes by this factor before plotting (default 1.0)
- `ylabel`: y-axis label; accepts any Makie-compatible string, including `LaTeXString`
"""
function plot_layer_volumes_timeseries(
    snaps::Vector{CO2BatchFill.LayerSnapshot};
    output_file::Union{String,Nothing}=nothing,
    vol_scale::Float64=1.0,
    ylabel="Volume",
    show_injected::Bool=false,
    linewidth::Int=_LW,
)
    title = isempty(snaps) ? "Layer Volumes" : snaps[1].layer_name
    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1]; xlabel="Time", ylabel, title)
    _timeseries_panel!(ax, snaps; show_injected, linewidth, vol_scale)
    axislegend(ax; patchsize=(40, 10), position=:lt)

    if output_file !== nothing
        save(output_file, fig)
        println("Saved layer volumes timeseries to: $(output_file)")
    end
    return fig
end

"""
Plot CO2 volume time-series for a multi-layer simulation.

- `vol_scale`: multiply raw SWIM volumes by this factor before plotting (default 1.0)
- `ylabel`: y-axis label; accepts any Makie-compatible string, including `LaTeXString`
"""
function plot_multi_layer_volumes_timeseries(
    snaps::Vector{CO2BatchFill.MultiLayerSnapshot};
    output_file::String="multi_layer_timeseries.svg",
    vol_scale::Float64=1.0,
    ylabel="Volume",
    show_injected::Bool=false,
    linewidth::Int=_LW,
    figure_size::Union{Tuple{Int,Int},Nothing}=nothing,
    figure_layout::Union{Tuple{Int,Int},Nothing}=nothing,
)
    n = length(snaps[1].layers)

    figure_size = isnothing(figure_size) ? (500 * n, 400) : figure_size
    figure_layout = isnothing(figure_layout) ? (1, n) : figure_layout

    fig = Figure(size=figure_size)
    axes = Axis[]
    for i in 1:n
        layer_snaps = [s.layers[i] for s in snaps]
        figure_row, figure_col = i_to_layout(figure_layout, i)
        ax = Axis(fig[figure_row, figure_col];
            xlabel="Time",
            ylabel=(i == 1 ? ylabel : ""),
            title=layer_snaps[1].layer_name,
        )
        _timeseries_panel!(ax, layer_snaps; show_injected, linewidth, vol_scale)
        if i == n
            axislegend(ax; patchsize=(40, 10), position=:lt)
        end
        push!(axes, ax)
    end
    linkyaxes!(axes...)

    if output_file !== nothing
        save(output_file, fig)
        println("Saved multi-layer volumes timeseries to: $(output_file)")
    end
    return fig
end

"""
Create an animated bird's eye view of CO2 filling in a single layer.
"""
function animate_layer_filling(
    layer::CO2BatchFill.Layer,
    seq::Vector{SpillEvent},
    domain::CO2BatchFill.Domain3D;
    output_file::String="layer_filling.gif",
    num_frames::Int=60,
    start_time::Union{Float64,Nothing}=nothing,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=15,
    colormap::Symbol=_CMAP,
    max_co2_height::Float64=_MAXH,
    show_contours::Bool=true,
    contour_levels::Int=10,
)
    tstruct = layer.trap_structure
    pad = layer.pad_width
    nx_padded, ny_padded = get_padded_grid_size(layers)
    nx, ny = get_grid_size(layers)

    t_start = isnothing(start_time) ? seq[1].timestamp : start_time
    if isnothing(end_time)
        finite_times = [se.timestamp for se in seq if isfinite(se.timestamp)]
        t_end = isempty(finite_times) ? t_start + 1.0 : maximum(finite_times)
    else
        t_end = end_time
    end

    timepoints = collect(range(t_start, stop=t_end, length=num_frames))

    println("Computing trap states for $(num_frames) frames...")
    tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)
    z_vol_tables = _compute_z_vol_tables(tstruct)

    topo = tstruct.topography[pad+1:end-pad, pad+1:end-pad]
    x_coords = range(0.0, nx * domain.dx, length=nx)
    y_coords = range(0.0, ny * domain.dy, length=ny)

    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1];
        xlabel="x (m)",
        ylabel="y (m)",
        title=layer.name,
        aspect=DataAspect(),
    )

    height_data = Observable(zeros(Float64, nx, ny))
    heatmap!(ax, x_coords, y_coords, height_data;
        colormap=colormap,
        colorrange=(0.0, max_co2_height),
    )
    Colorbar(fig[1, 2]; colormap, colorrange=(0.0, max_co2_height), label="Column height")

    if show_contours
        contour!(ax, x_coords, y_coords, topo;
            levels=contour_levels,
            color=(:white),
            linewidth=0.8,
        )
    end


    # Running-max volumes: prevent drainage from causing colors to decrease
    running_max_vols = zeros(Float64, numtraps(tstruct))

    println("Recording animation to: $(output_file)")
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]
        _, volumes, _ = tstates[frame_idx]
        running_max_vols .= max.(running_max_vols, volumes)
        height_data[] = _height_map(tstruct, z_vol_tables, running_max_vols, pad, nx_padded, ny_padded)
        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  frame $(frame_idx)/$(length(timepoints))")
        end
    end


    return fig
end

"""
Create an animated bird's eye view of CO2 filling across multiple layers,
with N side-by-side panels and a shared time label.
"""
function animate_multi_layer_filling(
    layers::Vector{CO2BatchFill.Layer},
    seqs::Vector{Vector{SpillEvent}},
    domain::CO2BatchFill.Domain3D;
    output_file::String="multi_layer_filling.gif",
    num_frames::Int=60,
    start_time::Union{Float64,Nothing}=nothing,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=15,
    colormap::Symbol=_CMAP,
    max_co2_height::Float64=_MAXH,
    show_contours::Bool=true,
    contour_levels::Int=10,
)
    n = length(layers)
    @assert n == length(seqs) "Must have one seq per layer"

    # Time range = union across all layers
    if isnothing(start_time)
        t_start = minimum(seq[1].timestamp for seq in seqs)
    else
        t_start = start_time
    end
    if isnothing(end_time)
        all_finite = [se.timestamp for seq in seqs for se in seq if isfinite(se.timestamp)]
        t_end = isempty(all_finite) ? t_start + 1.0 : maximum(all_finite)
    else
        t_end = end_time
    end

    timepoints = collect(range(t_start, stop=t_end, length=num_frames))

    println("Computing trap states for $(num_frames) frames per layer...")
    layer_data = map(1:n) do i
        tstruct = layers[i].trap_structure
        pad = layers[i].pad_width
        nx_p, ny_p = size(tstruct.topography)
        tstates = trap_states_at_timepoints(tstruct, seqs[i], timepoints; verbose=false)
        z_vol_tables = _compute_z_vol_tables(tstruct)
        (tstruct=tstruct, pad=pad, nx_p=nx_p, ny_p=ny_p,
            tstates=tstates, z_vol_tables=z_vol_tables)
    end

    fig = Figure(size=(700 * n, 650))

    height_observables = Vector{Observable{Matrix{Float64}}}(undef, n)
    for i in 1:n
        ld = layer_data[i]
        nx = ld.nx_p - 2 * ld.pad
        ny = ld.ny_p - 2 * ld.pad
        x_coords = range(0.0, nx * domain.dx, length=nx)
        y_coords = range(0.0, ny * domain.dy, length=ny)
        topo = ld.tstruct.topography[ld.pad+1:end-ld.pad, ld.pad+1:end-ld.pad]

        ax = Axis(fig[1, i];
            xlabel="x (m)",
            ylabel=(i == 1 ? "y (m)" : ""),
            title=layers[i].name,
            aspect=DataAspect(),
        )
        height_observables[i] = Observable(zeros(Float64, nx, ny))
        heatmap!(ax, x_coords, y_coords, height_observables[i];
            colormap=colormap, colorrange=(0.0, max_co2_height))
        if show_contours
            contour!(ax, x_coords, y_coords, topo;
                levels=contour_levels, color=(:white, 0.5), linewidth=0.8)
        end
    end
    Colorbar(fig[1, n+1]; colormap, colorrange=(0.0, max_co2_height), label="Column height")

    # Running-max volumes per layer: prevent drainage from causing colors to decrease
    running_max_vols = [zeros(Float64, numtraps(layer_data[i].tstruct)) for i in 1:n]

    println("Recording animation to: $(output_file)")
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]
        for i in 1:n
            ld = layer_data[i]
            _, volumes, _ = ld.tstates[frame_idx]
            running_max_vols[i] .= max.(running_max_vols[i], volumes)
            height_observables[i][] = _height_map(
                ld.tstruct, ld.z_vol_tables, running_max_vols[i], ld.pad, ld.nx_p, ld.ny_p)
        end
        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  frame $(frame_idx)/$(length(timepoints))")
        end
    end

    return fig
end

"""
Draw leakage location markers (triangles) for a single layer panel.
"""
function _draw_leakage_locations!(ax, layer, layer_snap, domain)
    tstruct = layer.trap_structure
    pad = layer.pad_width
    leak_xs, leak_ys = Float64[], Float64[]
    for trap_id in 1:numtraps(tstruct)
        layer_snap.trap_leaking[trap_id] || continue
        ci = find_leakage_location(trap_id, tstruct)
        x = (ci[1] - 1 - pad) * domain.dx
        y = (ci[2] - 1 - pad) * domain.dy
        push!(leak_xs, x)
        push!(leak_ys, y)
    end
    if !isempty(leak_xs)
        scatter!(ax, leak_xs, leak_ys;
            color=:black, marker=:utriangle, markersize=30)
    end
end

"""
Plot CO2 plume outlines for an ensemble of multi-layer simulations, with one
coloured contour line per ensemble member per layer panel.
"""
function plot_multi_layer_ensemble(
    layers::Vector{CO2BatchFill.Layer},
    ensemble::Vector{CO2BatchFill.MultiLayerSnapshot},
    domain::CO2BatchFill.Domain3D;
    contour_color=nothing,
    output_file::Union{String,Nothing}=nothing,
    contour_level::Float64=0.01,
    linewidth::Real=2.5,
    show_topography::Bool=true,
    show_labels::Bool=false,
    topo_contour_levels::Int=10,
    major_contour_every::Int=5,
    contour_opacity::Float64=0.8,
    injection_locations::Union{Nothing,Vector{Tuple{Float64,Float64}}}=nothing,
    show_leakage_locations::Bool=false,
    figure_size::Union{Tuple{Int,Int},Nothing}=nothing,
    member_labels::Union{Nothing,Vector{String}}=nothing,
    member_label_fontsize::Real=12,
)
    n_layers = length(layers)

    eff_color = isnothing(contour_color) ? :black : contour_color

    figure_size = isnothing(figure_size) ? (700 * n_layers, 600) : figure_size
    fig = Figure(size=figure_size)

    for i in 1:n_layers
        layer = layers[i]
        tstruct = layer.trap_structure
        pad = layer.pad_width
        nx_padded, ny_padded = get_padded_grid_size(layers)
        nx, ny = get_grid_size(layers)

        x_coords = range(0.0, nx * domain.dx, length=nx)
        y_coords = range(0.0, ny * domain.dy, length=ny)
        topo = tstruct.topography[pad+1:end-pad, pad+1:end-pad]
        z_vol_tables = _compute_z_vol_tables(tstruct)

        layer_name = ensemble[1].layers[i].layer_name
        ax = Axis(fig[1, i];
            xlabel="x",
            ylabel=(i == 1 ? "y" : ""),
            title=layer_name,
            aspect=DataAspect(),
        )

        if show_topography
            _topo_contours!(ax, x_coords, y_coords, topo;
                contour_levels=topo_contour_levels, major_contour_every,
                contour_opacity, show_labels)
        end

        for (j, snap) in enumerate(ensemble)
            layer_snap = snap.layers[i]
            height_array = _height_map(
                tstruct, z_vol_tables, layer_snap.trap_volumes, pad, nx_padded, ny_padded;
                leaking=layer_snap.trap_leaking,
                leakage_heights=layer_snap.trap_leakage_height)

            any(>=(contour_level), height_array) || continue

            has_label = !isnothing(member_labels)
            label_text = has_label ? member_labels[j] : ""
            contour!(ax, x_coords, y_coords, height_array;
                levels=[contour_level],
                color=eff_color,
                linewidth=linewidth,
                labels=has_label,
                labelsize=member_label_fontsize,
                labelformatter=(val -> label_text),
            )
        end

        if !isnothing(injection_locations) && i == 1
            scatter!(ax,
                [loc[1] for loc in injection_locations],
                [loc[2] for loc in injection_locations];
                color=:black, marker=:xcross, markersize=35)
        end

        if show_leakage_locations
            _draw_leakage_locations!(ax, layer, ensemble[end].layers[i], domain)
        end
    end

    if output_file !== nothing
        save(output_file, fig)
        println("Saved ensemble plot to: $(output_file)")
    end

    return fig
end

"""
Plot CO2 volume time-series with ensemble mean and uncertainty bands for a
multi-layer simulation.
"""
function plot_multi_layer_ensemble_timeseries(
    ensemble::Vector{Vector{CO2BatchFill.MultiLayerSnapshot}};
    output_file::Union{String,Nothing}=nothing,
    vol_scale::Float64=1.0,
    ylabel="Volume",
    linewidth::Int=_LW,
    figure_size::Union{Tuple{Int,Int},Nothing}=nothing,
)
    n_layers = length(ensemble[1][1].layers)
    n_ensemble = length(ensemble)
    n_times = length(ensemble[1])

    figure_size = isnothing(figure_size) ? (500 * n_layers, 400) : figure_size
    fig = Figure(size=figure_size)
    axes = Axis[]

    times = [s.timestamp for s in ensemble[1]]

    for i in 1:n_layers
        ax = Axis(fig[1, i];
            xlabel="Time",
            ylabel=(i == 1 ? ylabel : ""),
            title=ensemble[1][1].layers[i].layer_name,
        )

        stored_matrix = zeros(n_ensemble, n_times)
        drained_matrix = zeros(n_ensemble, n_times)
        for (j, member) in enumerate(ensemble)
            for (k, snap) in enumerate(member)
                stored_matrix[j, k] = snap.layers[i].total_stored * vol_scale
                drained_matrix[j, k] = snap.layers[i].total_drained * vol_scale
            end
        end

        stored_mean = vec(mean(stored_matrix; dims=1))
        stored_std = vec(std(stored_matrix; dims=1))
        drained_mean = vec(mean(drained_matrix; dims=1))
        drained_std = vec(std(drained_matrix; dims=1))

        band!(ax, times, stored_mean .- stored_std, stored_mean .+ stored_std;
            color=(:black, 0.15))
        lines!(ax, times, stored_mean;
            label="Stored", color=:black, linewidth=linewidth, linestyle=:solid)

        band!(ax, times, drained_mean .- drained_std, drained_mean .+ drained_std;
            color=(:black, 0.15))
        lines!(ax, times, drained_mean;
            label="Drained", color=:black, linewidth=linewidth, linestyle=:dash)

        if i == n_layers
            axislegend(ax; patchsize=(40, 10), position=:lt)
        end
        push!(axes, ax)
    end
    linkyaxes!(axes...)

    if output_file !== nothing
        save(output_file, fig)
        println("Saved ensemble plot to: $(output_file)")
    end

    return fig
end

end # module
