# CO2InjectionModeling

![Animation of CO2 Injection](media/animation.gif)

# CO2InjectionModeling.jl

CO2InjectionModeling.jl is a Julia package for modeling CO₂ injection and migration in subsurface geological formations. It provides tools to simulate the movement, trapping, and potential leakage of injected CO₂, with a focus on realistic geological structures and caprock topography.

## Features
- Layered geological models with customizable topography
- Simulation of CO2 migration, trapping, and leakage
- Visualization of injection scenarios and results
- Animation support for time-dependent processes

## Integration with SurfaceWaterIntegratedModeling
This package leverages the [SurfaceWaterIntegratedModeling.jl](https://github.com/your-org/SurfaceWaterIntegratedModeling.jl) package for advanced hydrological and surface water modeling. Key functions such as trap structure analysis, spill event handling, and time series interpolation are provided by SurfaceWaterIntegratedModeling and are used to:
- Identify and analyze geological traps
- Track CO₂ fill states and spill events over time
- Interpolate and visualize time-dependent results

## Getting Started
1. Clone the repository and install dependencies:
	```julia
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	```
2. Run example simulations from the `examples/` folder, e.g.:
	```julia
	include("examples/sleipner_layers.jl")
	```
3. View results and animations in the `media/` folder.

## Example Output
The package produces visualizations and animations, such as the one shown above, to illustrate CO2 migration and trapping over time.

## License
See [LICENSE](LICENSE).
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ellingsvee.github.io/CO2InjectionModeling.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ellingsvee.github.io/CO2InjectionModeling.jl/dev/)
[![Build Status](https://github.com/ellingsvee/CO2InjectionModeling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ellingsvee/CO2InjectionModeling.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ellingsvee/CO2InjectionModeling.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ellingsvee/CO2InjectionModeling.jl)
