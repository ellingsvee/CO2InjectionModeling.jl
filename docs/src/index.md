# CO2BatchFill.jl

CO2BatchFill.jl is an open-source Julia package for simulating horizontal and
vertical CO2 migration in stacked geological reservoirs.  It builds on the
[SurfaceWaterIntegratedModeling (SWIM)](https://github.com/sintefmath/SurfaceWaterIntegratedModeling.jl)
spill-analysis framework to provide:

- **Trap analysis** — identify structural traps from depth surfaces using SWIM's
  spill-point algorithm.
- **Single-layer filling** — simulate CO2 injection, structural trapping, and
  caprock leakage in one sand layer.
- **Multi-layer filling** — propagate CO2 leakage upward through a stack of
  layers, with per-layer injection, leakage detection, and residual trapping.
- **Visualization** — spatial CO2 height maps and volume time-series plots
  (requires a [Makie](https://makie.juliaplots.org/) backend).

## Contents

```@contents
Pages = ["types.md", "topography.md", "layer_analysis.md",
         "simulation.md", "analysis.md",
         "unit_conversion.md", "visualization.md", "utils.md"]
Depth = 2
```
