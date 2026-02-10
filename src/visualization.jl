# Visualization functions - implementations provided by MakieExt extension
# These functions are only available when Makie (or a Makie backend) is loaded

export animate_multi_layer_filling, plot_layer_volumes_timeseries
export plot_layer_topographies, plot_final_co2_distribution

"""
    animate_multi_layer_filling(layers, seqs, domain; kwargs...)

Create an animated bird's eye view of CO2 filling in all layers simultaneously.
Requires Makie (or a backend like CairoMakie) to be loaded.

See VISUALIZATION.md for full documentation and usage examples.
"""
function animate_multi_layer_filling end

"""
    plot_layer_volumes_timeseries(snapshots; kwargs...)

Plot CO2 volumes in each layer as a function of time using subplots.
Requires Makie (or a backend like CairoMakie) to be loaded.

See VISUALIZATION.md for full documentation and usage examples.
"""
function plot_layer_volumes_timeseries end

"""
    plot_layer_topographies(layers, domain; kwargs...)

Create a grid of contour plots with heatmaps showing the topography (depth) of each layer.
Requires Makie (or a backend like CairoMakie) to be loaded.

See VISUALIZATION.md for full documentation and usage examples.
"""
function plot_layer_topographies end

"""
    plot_final_co2_distribution(layers, seqs, domain; kwargs...)

Create a poster-ready static plot showing CO2 distribution at the final (or specified) timepoint.
Requires Makie (or a backend like CairoMakie) to be loaded.

See VISUALIZATION.md for full documentation and usage examples.
"""
function plot_final_co2_distribution end
