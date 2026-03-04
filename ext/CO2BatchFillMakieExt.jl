module CO2BatchFillMakieExt

using CO2BatchFill
using Makie
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints,
    SpillEvent, _compute_z_vol_tables

import CO2BatchFill: animate_layer_filling, animate_multi_layer_filling,
    plot_layer, plot_multi_layer,
    plot_layer_volumes_timeseries, plot_multi_layer_volumes_timeseries

# ─────────────────────────────────────────────────────────────────────────────
# Shared style constants
# ─────────────────────────────────────────────────────────────────────────────

const _CMAP = :thermal
const _MAXH = 20.0    # default max CO₂ height (m)
const _LW = 2       # default linewidth
const _COLORS = (injected=:black, stored=:blue, drained=:orange, passthru=:green, to_next=:red)

# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
Build an unpadded 2-D CO₂ height array from per-trap volumes.
"""
function _height_map(tstruct, z_vol_tables, volumes, pad, nx_p, ny_p)
    num_traps = numtraps(tstruct)
    height_map_padded = zeros(Float64, nx_p, ny_p)
    for trap_id in 1:num_traps
        vol = volumes[trap_id]
        vol <= 0.0 && continue
        h = volume_to_height(vol, trap_id, z_vol_tables[trap_id], tstruct)
        h <= 0.0 && continue
        for idx in tstruct.footprints[trap_id]
            height_map_padded[idx] = max(height_map_padded[idx], h)
        end
    end
    return height_map_padded[pad+1:end-pad, pad+1:end-pad]
end

"""
Draw a CO₂ height heatmap + topography contours for one `LayerSnapshot` onto `ax`.
Returns the heatmap object.
"""
function _layer_panel!(ax, layer, snap, domain;
    pad_width=2, colormap=_CMAP, max_co2_height=_MAXH,
    show_contours=true, contour_levels=10)

    tstruct = layer.trap_structure
    pad = layer.boundary_condition == :closed ? pad_width : 0
    nx_padded, ny_padded = size(tstruct.topography)
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    z_vol_tables = _compute_z_vol_tables(tstruct)
    height_array = _height_map(tstruct, z_vol_tables, snap.trap_volumes, pad, nx_padded, ny_padded)

    x_coords = range(0.0, nx * domain.dx, length=nx)
    y_coords = range(0.0, ny * domain.dy, length=ny)
    topo = tstruct.topography[pad+1:end-pad, pad+1:end-pad]

    hm = heatmap!(ax, x_coords, y_coords, height_array;
        colormap=colormap,
        colorrange=(0.0, max_co2_height),
    )

    if show_contours
        contour!(ax, x_coords, y_coords, topo;
            levels=contour_levels,
            color=(:white, 0.5),
            linewidth=0.8,
        )
    end

    return hm
end

"""
Draw stored/drained/passthrough/→next-layer volume lines for a vector of
`LayerSnapshot`s onto `ax`.
"""
function _timeseries_panel!(ax, snaps::Vector{LayerSnapshot};
    show_injected=false, linewidth=_LW)

    times = [s.timestamp for s in snaps]

    if show_injected
        lines!(ax, times, [s.total_injected for s in snaps];
            label="Injected", color=_COLORS.injected, linewidth=linewidth)
    end
    lines!(ax, times, [s.total_stored for s in snaps];
        label="Stored", color=_COLORS.stored, linewidth=linewidth)
    lines!(ax, times, [s.total_drained for s in snaps];
        label="Drained", color=_COLORS.drained, linewidth=linewidth)
    # lines!(ax, times, [s.total_passthrough for s in snaps];
    #     label="Passthrough", color=_COLORS.passthru, linewidth=linewidth)
    # lines!(ax, times, [total_to_next_layer(s) for s in snaps];
    #     label="→ Next layer", color=_COLORS.to_next, linewidth=linewidth)
end

# ─────────────────────────────────────────────────────────────────────────────
# Static spatial plots
# ─────────────────────────────────────────────────────────────────────────────

"""
Plot CO₂ height distribution for a single layer at a given snapshot time.
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
        xlabel="x (m)",
        ylabel="y (m)",
        title=snap.layer_name,
        aspect=DataAspect(),
    )
    _layer_panel!(ax, layer, snap, domain;
        pad_width, colormap, max_co2_height, show_contours, contour_levels)
    Colorbar(fig[1, 2]; colormap, colorrange=(0.0, max_co2_height), label="CO₂ height (m)")

    save(output_file, fig)
    println("Saved layer plot to: $(output_file)")
    return nothing
end

"""
Plot CO₂ height distribution for multiple layers side-by-side with a shared colorbar.
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
    contour_levels::Int=10,
)
    n = length(layers)
    fig = Figure(size=(700 * n, 600))

    for i in 1:n
        ax = Axis(fig[1, i];
            xlabel="x (m)",
            ylabel=(i == 1 ? "y (m)" : ""),
            title=snap.layers[i].layer_name,
            aspect=DataAspect(),
        )
        _layer_panel!(ax, layers[i], snap.layers[i], domain;
            pad_width, colormap, max_co2_height, show_contours, contour_levels)
    end
    Colorbar(fig[1, n+1]; colormap, colorrange=(0.0, max_co2_height), label="CO₂ height (m)")

    save(output_file, fig)
    println("Saved multi-layer plot to: $(output_file)")
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Time-series plots
# ─────────────────────────────────────────────────────────────────────────────

"""
Plot CO₂ volume time-series for a single layer.
"""
function plot_layer_volumes_timeseries(
    snaps::Vector{CO2BatchFill.LayerSnapshot};
    output_file::String="layer_timeseries.svg",
    show_injected::Bool=true,
    linewidth::Int=_LW,
)
    title = isempty(snaps) ? "Layer Volumes" : snaps[1].layer_name
    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1];
        xlabel="Time",
        ylabel="Volume (SWIM units)",
        title=title,
    )
    _timeseries_panel!(ax, snaps; show_injected, linewidth)
    axislegend(ax)

    save(output_file, fig)
    println("Saved layer volumes timeseries to: $(output_file)")
    return nothing
end

"""
Plot CO₂ volume time-series for a multi-layer simulation.

- `mode=:per_layer` (default): one panel per layer side-by-side; legend on last panel
- `mode=:system`: single panel with system totals (Injected / Stored / Surface leakage)
"""
function plot_multi_layer_volumes_timeseries(
    snaps::Vector{CO2BatchFill.MultiLayerSnapshot};
    output_file::String="multi_layer_timeseries.svg",
    mode::Symbol=:per_layer,
    show_injected::Bool=true,
    linewidth::Int=_LW,
)
    if mode == :per_layer
        n = length(snaps[1].layers)
        fig = Figure(size=(500 * n, 400))
        for i in 1:n
            layer_snaps = [s.layers[i] for s in snaps]
            ax = Axis(fig[1, i];
                xlabel="Time",
                ylabel=(i == 1 ? "Volume (SWIM units)" : ""),
                title=layer_snaps[1].layer_name,
            )
            _timeseries_panel!(ax, layer_snaps; show_injected, linewidth)
            if i == n
                axislegend(ax)
            end
        end
    elseif mode == :system
        fig = Figure(size=(700, 500))
        ax = Axis(fig[1, 1];
            xlabel="Time",
            ylabel="Volume (SWIM units)",
            title="System Totals",
        )
        times = [s.timestamp for s in snaps]
        if show_injected
            lines!(ax, times, [s.total_injected for s in snaps];
                label="Injected", color=_COLORS.injected, linewidth=linewidth)
        end
        lines!(ax, times, [s.total_stored for s in snaps];
            label="Stored", color=_COLORS.stored, linewidth=linewidth)
        lines!(ax, times, [s.total_surface_leakage for s in snaps];
            label="Surface leakage", color=_COLORS.to_next, linewidth=linewidth)
        axislegend(ax)
    else
        error("Unknown mode: $(mode). Use :per_layer or :system.")
    end

    save(output_file, fig)
    println("Saved multi-layer volumes timeseries to: $(output_file)")
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Animations
# ─────────────────────────────────────────────────────────────────────────────

"""
Create an animated bird's eye view of CO₂ filling in a single layer.
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
    Colorbar(fig[1, 2]; colormap, colorrange=(0.0, max_co2_height), label="CO₂ height (m)")

    if show_contours
        contour!(ax, x_coords, y_coords, topo;
            levels=contour_levels,
            color=(:white, 0.5),
            linewidth=0.8,
        )
    end

    time_label = Observable("t = $(round(t_start, digits=2))")
    Label(fig[0, :], time_label; fontsize=14)

    println("Recording animation to: $(output_file)")
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]
        _, volumes, _ = tstates[frame_idx]
        height_data[] = _height_map(tstruct, z_vol_tables, volumes, pad, nx_padded, ny_padded)
        time_label[] = "t = $(round(tp, digits=2))"
        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Done.")
    return nothing
end

"""
Create an animated bird's eye view of CO₂ filling across multiple layers,
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
    time_label = Observable("t = $(round(t_start, digits=2))")
    Label(fig[0, 1:n], time_label; fontsize=14)

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
    Colorbar(fig[1, n+1]; colormap, colorrange=(0.0, max_co2_height), label="CO₂ height (m)")

    println("Recording animation to: $(output_file)")
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]
        for i in 1:n
            ld = layer_data[i]
            _, volumes, _ = ld.tstates[frame_idx]
            height_observables[i][] = _height_map(
                ld.tstruct, ld.z_vol_tables, volumes, ld.pad, ld.nx_p, ld.ny_p)
        end
        time_label[] = "t = $(round(tp, digits=2))"
        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Done.")
    return nothing
end

end # module
