# Layer Analysis

Before running a simulation, each sand layer must be analysed to identify its
structural traps using the SWIM spill-point algorithm. [`analyze_base_surfaces`](@ref) iterates over the sand layers in the topography, optionally adds closed boundary walls, and calls
SWIM's `spillanalysis` for each surface.

```@docs
CO2BatchFill.analyze_base_surfaces
CO2BatchFill.add_boundary_wall
```
