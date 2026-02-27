# CO2BatchFill.jl

Efficient modeling of multi-layer CO2 migration in geological reservoirs. A fast, simple and interpretable alternative to full-physics models. Uses trap-based filling via [SurfaceWaterIntegratedModeling.jl](https://github.com/sintefmath/SurfaceWaterIntegratedModeling.jl) with caprock leakage and multi-layer coupling.

Further described in this [poster](https://ellingsvee.github.io/Files/Posters/GeiloWinterScool.pdf).

![CO2 distribution](./media/co2_distribution_opt2wells.svg)

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

## Quick Start

A full working example is provided in [`examples/synthetic_example.jl`](examples/synthetic_example.jl). The basic workflow is:

```julia
using CO2BatchFill
using SurfaceWaterIntegratedModeling

# 1. Define topography (depth surfaces for each sand layer)
#    sand_layers ordered shallowest-first, deepest-last
sand_layers = [
    Dict{String,Any}("name" => "L2", "top" => surface_shallow, "base" => base_shallow),
    Dict{String,Any}("name" => "L1", "top" => surface_deep,    "base" => base_deep),
]
topography = GenericTopography(sand_layers, nx, ny, dx, dy, depth_min, depth_max)

# 2. Analyze trap structures
domain = create_domain(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=:closed)

# 3. Set reservoir properties
rprops = ReservoirProperties(
    0.30,      # porosity
    0.20,      # residual CO2 saturation
    0.10,      # irreducible water saturation
    20000.0,   # shale pressure threshold (Pa)
    5.0,       # residual leakage time (years)
)

# 4. Define injection (one event vector per layer)
injection_events = [layer1_events, layer2_events]

# 5. Run simulation
seqs, leakage_states = fill_layers(layers, domain, rprops, injection_events)

# 6. Analyze results
snapshots = generate_reservoir_snapshots(
    layers, seqs, leakage_states, domain, rprops, injection_events;
    num_snapshots=10, start_time=0.0, end_time=15.0,
)
print_snapshot_summary(snapshots[end])
```

## Testing

```julia
using Pkg
Pkg.test("CO2BatchFill")
```

## Visualization

Plotting is available via a Makie extension. Load any Makie backend to enable it:

```julia
using CairoMakie
plot_final_co2_distribution(layers, seqs, domain; time=15.0)
```
