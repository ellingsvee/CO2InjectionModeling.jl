# Topography

CO2BatchFill uses a topography abstraction to describe the geometry of the
reservoir.  The [`GenericTopography`](@ref) concrete type is suitable for most
use cases; custom types can be created by subtyping [`AbstractTopography`](@ref)
and implementing the accessor interface.

## Abstract interface

```@docs
CO2BatchFill.AbstractTopography
CO2BatchFill.get_sand_layers
CO2BatchFill.get_grid_dimensions
CO2BatchFill.get_grid_spacing
CO2BatchFill.get_depth_range
CO2BatchFill.get_caprock_surface
CO2BatchFill.get_coordinate_origin
CO2BatchFill.get_num_layers
```

## Concrete type

```@docs
CO2BatchFill.GenericTopography
```

## Domain creation

```@docs
CO2BatchFill.create_domain
```
