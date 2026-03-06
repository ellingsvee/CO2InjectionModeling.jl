# Visualization

Visualization functions are provided as a [Makie](https://makie.juliaplots.org/)
extension and are available only when a Makie backend has been loaded, e.g.:

```julia
using CairoMakie   # or GLMakie, WGLMakie, …
```

## Spatial plots

```@docs
CO2BatchFill.plot_multi_layer
CO2BatchFill.plot_layer
```

## Time-series plots

```@docs
CO2BatchFill.plot_multi_layer_volumes_timeseries
CO2BatchFill.plot_layer_volumes_timeseries
```

## Animations

```@docs
CO2BatchFill.animate_multi_layer_filling
CO2BatchFill.animate_layer_filling
```
