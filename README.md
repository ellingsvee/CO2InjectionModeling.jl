# CO2BatchFill.jl

Efficient modeling of multi-layer CO2 migration in geological reservoirs. A fast, simple and interpretable alternative to full-physics models. Uses trap-based filling via [SurfaceWaterIntegratedModeling.jl](https://github.com/sintefmath/SurfaceWaterIntegratedModeling.jl) with caprock leakage and multi-layer coupling.

Further described in this [poster](https://ellingsvee.github.io/Files/Posters/GeiloWinterScool.pdf).

Documentation is available at https://ellingsvee.github.io/CO2BatchFill.jl

## Installation

Requires Julia 1.10+. From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/ellingsvee/CO2BatchFill.jl")
```

Or for development:

```julia
Pkg.develop(url="https://github.com/ellingsvee/CO2BatchFill.jl")
```

## Usage

A full working example is provided in [`examples/multi_layer_filling.jl`](examples/multi_layer_filling.jl).

## Testing

```julia
using Pkg
Pkg.test("CO2BatchFill")
```

## Visualization

Plotting is available via a Makie extension. Load any Makie backend to enable it.
