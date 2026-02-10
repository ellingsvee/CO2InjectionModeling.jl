"""
    Sleipner

Example module for Sleipner CO2 storage site data and configuration.
Provides Sleipner-specific topography loading, reservoir properties,
and R interface convenience wrappers.

Usage:
```julia
include("examples/Sleipner/Sleipner.jl")
using .Sleipner
```
"""
module Sleipner

using CO2BatchFill
using CO2BatchFill: AbstractTopography, Domain3D, Layer, ReservoirProperties, InjectionEvent

include("topography.jl")
include("setup.jl")

end # module Sleipner
