# Visualization functions - implementations provided by MakieExt extension
# These functions are only available when Makie (or a Makie backend) is loaded

export animate_layer_filling, animate_multi_layer_filling
export plot_layer, plot_multi_layer, plot_multi_layer_ensemble
export plot_layer_volumes_timeseries, plot_multi_layer_volumes_timeseries,
    plot_multi_layer_ensemble_timeseries

"""
    animate_layer_filling(layer, seq, leakage_state, weather_events, timepoints, domain; kwargs...)

Create an animation of CO2 filling for a single layer over `timepoints`.

Requires a Makie backend (e.g. `using CairoMakie`).
"""
function animate_layer_filling end

"""
    animate_multi_layer_filling(layers, seqs, leakage_states, weather_events_per_layer,
                                timepoints, domain; kwargs...)

Create an animation of multi-layer CO2 filling over `timepoints`.

Requires a Makie backend (e.g. `using CairoMakie`).
"""
function animate_multi_layer_filling end

"""
    plot_layer(layer, snap, domain; kwargs...) -> nothing

Plot the CO2 column-height distribution for a single layer at the time captured
in `snap` and save to a file.

Requires a Makie backend (e.g. `using CairoMakie`).

# Arguments
- `layer`: [`Layer`](@ref) from [`analyze_base_surfaces`](@ref)
- `snap`: [`LayerSnapshot`](@ref) from [`generate_layer_snapshot`](@ref)
- `domain`: [`Domain3D`](@ref)

# Keyword arguments
- `output_file`: Path for the saved figure (default `"layer_co2.svg"`)
- `pad_width`: Boundary padding cells (default `2`)
- `colormap`: Makie colormap (default `:thermal`)
- `max_co2_height`: Colorbar upper limit in metres (default `20.0`)
- `show_contours`: Overlay topography contours (default `true`)
- `contour_levels`: Number of contour levels or explicit vector (default `10`)
"""
function plot_layer end

"""
    plot_multi_layer(layers, snap, domain; kwargs...) -> nothing

Plot CO2 column-height distributions for all layers side-by-side with a shared
colorbar and save to a file.

Requires a Makie backend (e.g. `using CairoMakie`).

# Arguments
- `layers`: `Vector{Layer}` from [`analyze_base_surfaces`](@ref)
- `snap`: [`MultiLayerSnapshot`](@ref)
- `domain`: [`Domain3D`](@ref)

# Keyword arguments
- `output_file`: Output file path (default `"multi_layer_co2.svg"`)
- `pad_width`: Boundary padding (default `2`)
- `colormap`: Makie colormap symbol (default `:thermal`)
- `max_co2_height`: Colorbar upper limit in metres (default `20.0`)
- `show_contours`: Overlay topography contours (default `true`)
- `show_labels`: Label major contours (default `false`)
- `contour_levels`: Contour count or vector (default `10`)
- `major_contour_every`: Every nth level is drawn as a major contour (default `5`)
- `contour_opacity`: Contour line opacity (default `0.8`)
- `injection_locations`: `Vector{Tuple{Float64,Float64}}` of (x, y) well locations
  to mark on layer 1 (default `nothing`)
- `show_leakage_locations`: Mark leakage points with triangles (default `false`)
- `figure_size`: `(width, height)` in pixels (default `(700*n_layers, 600)`)
"""
function plot_multi_layer end

"""
    plot_multi_layer_ensemble(layers, ensemble, domain; kwargs...) -> nothing

Plot CO2 plume outlines for an ensemble of multi-layer simulations.
Each ensemble member is drawn as a single contour line (at `contour_level` metres of
CO2 column height) in a distinct colour, giving a compact uncertainty overview.
Panels are arranged side-by-side, one per geological layer.

Requires a Makie backend (e.g. `using CairoMakie`).

# Arguments
- `layers`: `Vector{Layer}` from [`analyze_base_surfaces`](@ref)
- `ensemble`: `Vector{MultiLayerSnapshot}`, one element per ensemble member
- `domain`: [`Domain3D`](@ref)

# Keyword arguments
- `output_file`: Output file path (default `"ensemble_co2.svg"`)
- `labels`: `Vector{String}` legend labels, one per member (default `"Member i"`)
- `colors`: Override colours (default: Wong colour palette, cycling if needed)
- `pad_width`: Boundary padding cells (default `2`)
- `contour_level`: CO2 column height threshold that defines each plume outline in metres
  (default `0.01`)
- `linewidth`: Outline line width (default `2.5`)
- `show_topography`: Overlay topography contours (default `true`)
- `topo_contour_levels`: Number of topography contour levels (default `10`)
- `major_contour_every`: Every nth level is drawn as a major contour (default `5`)
- `injection_locations`: `Vector{Tuple{Float64,Float64}}` of (x, y) well locations
  to mark on layer 1 (default `nothing`)
- `figure_size`: `(width, height)` in pixels (default `(700*n_layers, 600)`)
"""
function plot_multi_layer_ensemble end

"""
    plot_layer_volumes_timeseries(snaps; kwargs...) -> nothing

Plot stored and drained CO2 volumes over time for a single layer and save to a
file.

Requires a Makie backend (e.g. `using CairoMakie`).

# Arguments
- `snaps`: `Vector{LayerSnapshot}` from [`generate_layer_snapshots`](@ref)

# Keyword arguments
- `output_file`: Output file path (default `"layer_timeseries.svg"`)
- `vol_scale`: Multiply volumes by this factor before plotting (default `1.0`)
- `ylabel`: Y-axis label string or `LaTeXString` (default `"Volume"`)
- `show_injected`: Also plot injected volume (default `false`)
- `linewidth`: Line width (default `4`)
"""
function plot_layer_volumes_timeseries end

"""
    plot_multi_layer_volumes_timeseries(snaps; kwargs...) -> nothing

Plot CO2 volume time-series for every layer in a multi-layer simulation
side-by-side and save to a file.

Requires a Makie backend (e.g. `using CairoMakie`).

# Arguments
- `snaps`: `Vector{MultiLayerSnapshot}` from [`generate_multi_layer_snapshots`](@ref)

# Keyword arguments
- `output_file`: Output file path (default `"multi_layer_timeseries.svg"`)
- `vol_scale`: Volume scaling factor (default `1.0`)
- `ylabel`: Y-axis label (default `"Volume"`)
- `show_injected`: Also plot injected volume (default `false`)
- `linewidth`: Line width (default `4`)
- `figure_size`: `(width, height)` in pixels (default `(500*n_layers, 400)`)
"""
function plot_multi_layer_volumes_timeseries end

"""
    plot_multi_layer_ensemble_timeseries(ensemble; kwargs...) -> nothing

Plot CO2 volume time-series with ensemble mean and ±1σ uncertainty bands for
every layer in a multi-layer simulation, side-by-side.

Requires a Makie backend (e.g. `using CairoMakie`).

# Arguments
- `ensemble`: `Vector{Vector{MultiLayerSnapshot}}`, outer index = ensemble member,
  inner index = timepoint (all members must share the same timepoints)

# Keyword arguments
- `output_file`: Output file path (default `"ensemble_timeseries.svg"`)
- `vol_scale`: Volume scaling factor (default `1.0`)
- `ylabel`: Y-axis label (default `"Volume"`)
- `linewidth`: Line width (default `4`)
- `figure_size`: `(width, height)` in pixels (default `(500*n_layers, 400)`)
"""
function plot_multi_layer_ensemble_timeseries end