# CO2BatchFill.jl

[![CI](https://github.com/ellingsvee/CO2BatchFill.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ellingsvee/CO2BatchFill.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://ellingsvee.github.io/CO2BatchFill.jl)

A Julia package for simulating CO2 migration and structural trapping in stacked geological reservoirs. `CO2BatchFill` provides a fast, interpretable alternative to full-physics reservoir simulators by combining trap-based spill analysis from [SurfaceWaterIntegratedModeling.jl](https://github.com/sintefmath/SurfaceWaterIntegratedModeling.jl) (SWIM) with multi-layer coupling.

<p align="center">
  <img src="media/co2_distribution_opt2wells.svg" width="700"/>
</p>

## Key Features

- **Structural trap analysis:** Wraps algorithms from SWIM for identification of spill points. Handles both open and closed boundary conditions.
- **Multi-layer filling simulation:** CO2 injected into a deep layer fills traps, and eventually migrates upward through a stack of reservoir layers.
- **Visualization:**  A [Makie](https://docs.makie.org/) extension enables plotting of CO2 plumes and timeseries directly from `CO2BatchFill` data structures.

## Installation

Requires Julia 1.11+. From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/ellingsvee/CO2BatchFill.jl")
```


## Examples
There are currently three examples available in the `examples/` directory, which are also included in the documentation:

- [`multi_layer_filling.jl`](examples/multi_layer_filling.jl): Core workflow using a synthetic 3-layer domain with injection, migration, and plotting
- [`multi_layer_ensemble.jl`](examples/multi_layer_ensemble.jl): 100-run Monte Carlo ensemble varying capillary entry pressure.
- [`sleipner_analysis.jl`](examples/sleipner_analysis.jl): Real-world analysis using Sleipner depth surfaces.

## Testing

```julia
using Pkg
Pkg.test("CO2BatchFill")
```

## Documentation

Full API documentation is available at [ellingsvee.github.io/CO2BatchFill.jl](https://ellingsvee.github.io/CO2BatchFill.jl).


## License and contributions
`CO2BatchFill.jl` is licensed under the MIT License (see [LICENSE](LICENSE)). Contributions are welcome via pull requests or issues on GitHub. See the [CONTRIBUTING](CONTRIBUTING.md) guide for more details.