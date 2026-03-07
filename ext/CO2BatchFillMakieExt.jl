module CO2BatchFillMakieExt

using CO2BatchFill
using Makie
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints,
    SpillEvent, _compute_z_vol_tables
using Statistics: mean

import CO2BatchFill: animate_layer_filling, animate_multi_layer_filling,
    plot_layer, plot_multi_layer, plot_multi_layer_ensemble,
    plot_layer_volumes_timeseries, plot_multi_layer_volumes_timeseries

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
        for idx in tstruct.footprints[trap_id]
            height_map_padded[idx] = max(height_map_padded[idx], h)
        end
    end
    return height_map_padded[pad+1:end-pad, pad+1:end-pad]
end

"""
Draw a CO2 height heatmap + topography contours for one `LayerSnapshot` onto `ax`.
Returns the heatmap object.
"""
function _layer_panel!(ax, layer, snap, domain;
    pad_width=2, colormap=_CMAP, max_co2_height=_MAXH,
    show_contours=true, show_labels=false, contour_levels=10, major_contour_every=5, contour_opacity=0.8)

    tstruct = layer.trap_structure
    pad = layer.boundary_condition == :closed ? pad_width : 0
    nx_padded, ny_padded = size(tstruct.topography)
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    z_vol_tables = _compute_z_vol_tables(tstruct)
    height_array = _height_map(tstruct, z_vol_tables, snap.trap_volumes, pad, nx_padded, ny_padded;
        leaking=snap.trap_leaking,
        leakage_heights=snap.trap_leakage_height)

    x_coords = range(0.0, nx * domain.dx, length=nx)
    y_coords = range(0.0, ny * domain.dy, length=ny)
    topo = tstruct.topography[pad+1:end-pad, pad+1:end-pad]

    hm = heatmap!(ax, x_coords, y_coords, height_array;
        colormap=colormap,
        colorrange=(0.0, max_co2_height),
    )

    if show_contours
        # Determine all contour levels
        if contour_levels isa Integer
            mn, mx = extrema(topo)
            all_levels = collect(range(mn, mx; length=contour_levels))
        else
            all_levels = collect(contour_levels)
        end

        # Split into major / minor
        major_levels = all_levels[2:major_contour_every:end]
        minor_levels = setdiff(all_levels, major_levels)

        # Minor contours (thin, no labels)
        contour!(ax, x_coords, y_coords, topo;
            levels=minor_levels,
            color=(:black, contour_opacity),
            linewidth=0.8,
            labels=false
        )

        # Major contours (thicker, optional labels)
        contour!(ax, x_coords, y_coords, topo;
            levels=major_levels,
            color=(:black, contour_opacity),
            linewidth=2.0,
            labels=show_labels,
            labelsize=12,
        )
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
    output_file::String="layer_co2.svg",
    pad_width::Int=2,
    colormap::Symbol=_CMAP,
    max_co2_height::Float64=_MAXH,
    show_contours::Bool=true,
    contour_levels::Int=10,
)
    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1];
        xlabel="x",
        ylabel="y",
        title=snap.layer_name,
        aspect=DataAspect(),
    )
    _layer_panel!(ax, layer, snap, domain;
        pad_width, colormap, max_co2_height, show_contours, contour_levels)
    Colorbar(fig[1, 2]; colormap, colorrange=(0.0, max_co2_height), label="Column height", size=_CBAR_SIZE)

    save(output_file, fig)
    println("Saved layer plot to: $(output_file)")
    return nothing
end

"""
Plot CO2 height distribution for multiple layers side-by-side with a shared colorbar.
"""
function plot_multi_layer(
    layers::Vector{CO2BatchFill.Layer},
    snap::CO2BatchFill.MultiLayerSnapshot,
    domain::CO2BatchFill.Domain3D;
    output_file::String="multi_layer_co2.svg",
    pad_width::Int=2,
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
)
    n = length(layers)

    figure_size = isnothing(figure_size) ? (700 * n, 600) : figure_size
    fig = Figure(size=figure_size)

    has_injection = !isnothing(injection_locations)
    has_any_leakage = false

    for i in 1:n
        ax = Axis(fig[1, i];
            xlabel="x",
            ylabel=(i == 1 ? "y" : ""),
            title=snap.layers[i].layer_name,
            aspect=DataAspect(),
        )
        _layer_panel!(ax, layers[i], snap.layers[i], domain;
            pad_width, colormap, max_co2_height, show_contours, show_labels, contour_levels, major_contour_every, contour_opacity)

        if has_injection && i == 1
            scatter!(ax, [loc[1] for loc in injection_locations], [loc[2] for loc in injection_locations];
                color=:black, marker=:xcross, markersize=35)
        end

        if show_leakage_locations
            tstruct = layers[i].trap_structure
            pad = layers[i].boundary_condition == :closed ? pad_width : 0
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
    Colorbar(fig[1, n+1]; colormap, colorrange=(0.0, max_co2_height), label="Column height", size=_CBAR_SIZE)

    save(output_file, fig)
    println("Saved multi-layer plot to: $(output_file)")
    return nothing
end

# Time-series plots
"""
Plot CO2 volume time-series for a single layer.

- `vol_scale`: multiply raw SWIM volumes by this factor before plotting (default 1.0)
- `ylabel`: y-axis label; accepts any Makie-compatible string, including `LaTeXString`
"""
function plot_layer_volumes_timeseries(
    snaps::Vector{CO2BatchFill.LayerSnapshot};
    output_file::String="layer_timeseries.svg",
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

    save(output_file, fig)
    println("Saved layer volumes timeseries to: $(output_file)")
    return nothing
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
)
    n = length(snaps[1].layers)

    figure_size = isnothing(figure_size) ? (500 * n, 400) : figure_size
    fig = Figure(size=figure_size)
    axes = Axis[]
    for i in 1:n
        layer_snaps = [s.layers[i] for s in snaps]
        ax = Axis(fig[1, i];
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

    save(output_file, fig)
    println("Saved multi-layer volumes timeseries to: $(output_file)")
    return nothing
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
    pad_width::Int=2,
)
    tstruct = layer.trap_structure
    pad = layer.boundary_condition == :closed ? pad_width : 0
    nx_padded, ny_padded = size(tstruct.topography)
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

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

    println("Done.")
    return nothing
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
    pad_width::Int=2,
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
        pad = layers[i].boundary_condition == :closed ? pad_width : 0
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

    println("Done.")
    return nothing
end

# Wong (2011) colour palette — 7 colours that are readable by colour-blind viewers
const _WONG_COLORS = [
    RGBf(0.902, 0.624, 0.000),   # orange
    RGBf(0.337, 0.706, 0.914),   # sky blue
    RGBf(0.000, 0.620, 0.451),   # green
    RGBf(0.941, 0.894, 0.259),   # yellow
    RGBf(0.000, 0.447, 0.698),   # blue
    RGBf(0.835, 0.369, 0.000),   # vermillion
    RGBf(0.800, 0.475, 0.655),   # pink
]

"""
Plot CO2 plume outlines for an ensemble of multi-layer simulations, with one
coloured contour line per ensemble member per layer panel.
"""
function plot_multi_layer_ensemble(
    layers::Vector{CO2BatchFill.Layer},
    ensemble::Vector{CO2BatchFill.MultiLayerSnapshot},
    domain::CO2BatchFill.Domain3D;
    labels::Union{Nothing,Vector{String}}=nothing,
    colors=nothing,
    output_file::String="ensemble_co2.svg",
    pad_width::Int=2,
    contour_level::Float64=0.01,
    linewidth::Real=2.5,
    show_topography::Bool=true,
    topo_contour_levels::Int=10,
    major_contour_every::Int=5,
    injection_locations::Union{Nothing,Vector{Tuple{Float64,Float64}}}=nothing,
    figure_size::Union{Tuple{Int,Int},Nothing}=nothing,
)
    n_layers = length(layers)
    n_ensemble = length(ensemble)

    eff_labels = isnothing(labels) ? ["Member $i" for i in 1:n_ensemble] : labels
    eff_colors = if isnothing(colors)
        [_WONG_COLORS[((i - 1) % length(_WONG_COLORS)) + 1] for i in 1:n_ensemble]
    else
        colors
    end

    figure_size = isnothing(figure_size) ? (700 * n_layers, 600) : figure_size
    fig = Figure(size=figure_size)

    for i in 1:n_layers
        layer = layers[i]
        tstruct = layer.trap_structure
        pad = layer.boundary_condition == :closed ? pad_width : 0
        nx_padded, ny_padded = size(tstruct.topography)
        nx = nx_padded - 2 * pad
        ny = ny_padded - 2 * pad

        x_coords = range(0.0, nx * domain.dx, length=nx)
        y_coords = range(0.0, ny * domain.dy, length=ny)
        topo = tstruct.topography[pad+1:end-pad, pad+1:end-pad]
        z_vol_tables = _compute_z_vol_tables(tstruct)

        layer_name = ensemble[1].layers[i].layer_name
        ax = Axis(fig[1, i];
            xlabel="x (m)",
            ylabel=(i == 1 ? "y (m)" : ""),
            title=layer_name,
            aspect=DataAspect(),
        )

        if show_topography
            mn, mx = extrema(topo)
            all_levels = collect(range(mn, mx; length=topo_contour_levels))
            major_levels = all_levels[2:major_contour_every:end]
            minor_levels = setdiff(all_levels, major_levels)
            contour!(ax, x_coords, y_coords, topo;
                levels=minor_levels, color=(:black, 0.20), linewidth=0.6)
            contour!(ax, x_coords, y_coords, topo;
                levels=major_levels, color=(:black, 0.40), linewidth=1.2)
        end

        for (j, snap) in enumerate(ensemble)
            layer_snap = snap.layers[i]
            height_array = _height_map(
                tstruct, z_vol_tables, layer_snap.trap_volumes, pad, nx_padded, ny_padded;
                leaking=layer_snap.trap_leaking,
                leakage_heights=layer_snap.trap_leakage_height)

            any(>=(contour_level), height_array) || continue
            contour!(ax, x_coords, y_coords, height_array;
                levels=[contour_level],
                color=eff_colors[j],
                linewidth=linewidth,
            )
        end

        if !isnothing(injection_locations) && i == 1
            scatter!(ax,
                [loc[1] for loc in injection_locations],
                [loc[2] for loc in injection_locations];
                color=:black, marker=:xcross, markersize=35)
        end
    end

    # Build legend from LineElements — works regardless of which contours were drawn
    legend_elements = [LineElement(color=eff_colors[j], linewidth=linewidth) for j in 1:n_ensemble]
    Legend(fig[1, n_layers+1], legend_elements, eff_labels; framevisible=true)

    save(output_file, fig)
    println("Saved ensemble plot to: $(output_file)")
    return nothing
end

end # module
