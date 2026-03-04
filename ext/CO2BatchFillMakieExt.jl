module CO2BatchFillMakieExt

using CO2BatchFill
using Makie
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints,
    SpillEvent, _compute_z_vol_tables

import CO2BatchFill: animate_layer_filling


"""
    animate_layer_filling(layer, seq, domain; kwargs...)

Create an animated bird's eye view of CO2 filling in a single layer.

# Arguments
- `layer`  : `Layer` struct (from `analyze_base_surfaces`)
- `seq`    : `Vector{SpillEvent}` from `fill_sequence_with_leakage`
- `domain` : `Domain3D` struct

# Keyword Arguments
- `output_file`    : Path to save animation (default: `"layer_filling.mp4"`)
- `num_frames`     : Number of animation frames (default: `60`)
- `start_time`     : Start time (default: `seq[1].timestamp`)
- `end_time`       : End time; auto-detected from `seq` if `nothing` (default: `nothing`)
- `fps`            : Frames per second (default: `15`)
- `colormap`       : Colormap for CO2 height (default: `:thermal`)
- `max_co2_height` : Upper bound for colorscale in metres (default: `20.0`)
- `show_contours`  : Overlay terrain contour lines (default: `true`)
- `contour_levels` : Number of contour lines (default: `10`)
"""
function animate_layer_filling(
    layer::CO2BatchFill.Layer,
    seq::Vector{SpillEvent},
    domain::CO2BatchFill.Domain3D;
    output_file::String="layer_filling.mp4",
    num_frames::Int=60,
    start_time::Union{Float64,Nothing}=nothing,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=15,
    colormap::Symbol=:thermal,
    max_co2_height::Float64=20.0,
    show_contours::Bool=true,
    contour_levels::Int=10,
)
    tstruct = layer.trap_structure
    num_traps = numtraps(tstruct)
    pad = layer.boundary_condition == :closed ? 1 : 0

    nx_padded, ny_padded = size(tstruct.topography)
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # ── Determine time range ──────────────────────────────────────────────────
    t_start = isnothing(start_time) ? seq[1].timestamp : start_time

    if isnothing(end_time)
        finite_times = [se.timestamp for se in seq if isfinite(se.timestamp)]
        t_end = isempty(finite_times) ? t_start + 1.0 : maximum(finite_times)
    else
        t_end = end_time
    end

    timepoints = collect(range(t_start, stop=t_end, length=num_frames))

    # ── Precompute trap states and z→volume tables ────────────────────────────
    println("Computing trap states for $(num_frames) frames...")
    tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)
    z_vol_tables = _compute_z_vol_tables(tstruct)

    # ── Unpadded topography for contour overlay ───────────────────────────────
    topo = tstruct.topography[pad+1:end-pad, pad+1:end-pad]

    x_coords = range(0.0, nx * domain.dx, length=nx)
    y_coords = range(0.0, ny * domain.dy, length=ny)

    # ── Set up figure ─────────────────────────────────────────────────────────
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

    # ── Record animation ──────────────────────────────────────────────────────
    println("Recording animation to: $(output_file)")
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]
        filled, volumes, _ = tstates[frame_idx]

        height_map_padded = zeros(Float64, nx_padded, ny_padded)

        for trap_id in 1:num_traps
            vol = volumes[trap_id]
            (vol <= 0.0 && !filled[trap_id]) && continue

            h = volume_to_height(vol, trap_id, z_vol_tables[trap_id], tstruct)
            h <= 0.0 && continue

            for idx in tstruct.footprints[trap_id]
                height_map_padded[idx] = max(height_map_padded[idx], h)
            end
        end

        height_data[] = height_map_padded[pad+1:end-pad, pad+1:end-pad]
        time_label[] = "t = $(round(tp, digits=2))"

        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Done.")
    return nothing
end


end # module
