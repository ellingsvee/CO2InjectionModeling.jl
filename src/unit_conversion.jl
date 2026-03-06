using SurfaceWaterIntegratedModeling

export swim_volume_to_physical_volume, physical_volume_to_swim_volume
export full_volume_to_rock_volume_scaling, unit_volume_to_physical_scaling
export convert_injection_event_to_weather_event

"""
    unit_volume_to_physical_scaling(domain) -> Float64

Return the factor that converts a SWIM "unit volume" (one vertical unit per
horizontal cell) to physical volume in m³.  Equals `dx * dy`.
"""
function unit_volume_to_physical_scaling(domain::Domain3D)
    return domain.dx * domain.dy
end

"""
    full_volume_to_rock_volume_scaling(reservoir_properties) -> Float64

Return the pore-volume fraction available for CO2 storage:
`porosity * (1 - irreducible_water_saturation)`.
"""
function full_volume_to_rock_volume_scaling(
    reservoir_properties::ReservoirProperties,
)::Float64
    return reservoir_properties.sand_porosity * (1.0 - reservoir_properties.sand_irreducible_water_saturation)
end


"""
    swim_volume_to_physical_volume(volume, reservoir_properties, domain)

Convert SWIM internal volume units to physical volume in m³.

SWIM volumes are dimensionless (height × cell count); this function applies
the pore-volume scaling (`porosity × (1 - Swi)`) and the horizontal cell area
(`dx × dy`) to obtain a physically meaningful volume.

# Arguments
- `volume`: Scalar or array in SWIM units
- `reservoir_properties`: [`ReservoirProperties`](@ref)
- `domain`: [`Domain3D`](@ref)

# Example
```julia
scale = swim_volume_to_physical_volume(1.0, rp, domain)  # m³ per SWIM unit
physical_vol = snap.total_stored * scale
```
"""
function swim_volume_to_physical_volume(
    volume::Union{AbstractArray{<:Real},Real},
    reservoir_properties::ReservoirProperties,
    domain::Domain3D
)
    scaling_factor = full_volume_to_rock_volume_scaling(reservoir_properties) * unit_volume_to_physical_scaling(domain)
    return volume .* scaling_factor
end

"""
    physical_volume_to_swim_volume(volume, reservoir_properties, domain)

Convert physical volume in m³ to SWIM internal volume units.
Inverse of [`swim_volume_to_physical_volume`](@ref).
"""
function physical_volume_to_swim_volume(
    volume::Union{AbstractArray{<:Real},Real},
    reservoir_properties::ReservoirProperties,
    domain::Domain3D
)
    scaling_factor = full_volume_to_rock_volume_scaling(reservoir_properties) * unit_volume_to_physical_scaling(domain)
    return volume ./ scaling_factor
end

"""
    convert_injection_event_to_weather_event(injection_events, reservoir_properties, domain)
        -> Vector{WeatherEvent}

Convert a vector of [`InjectionEvent`](@ref)s (physical m³/yr rates) into SWIM
`WeatherEvent`s (dimensionless "rain rate" per cell) using the unit-conversion
scaling from [`physical_volume_to_swim_volume`](@ref).

Called internally by [`fill_layers`](@ref).
"""
function convert_injection_event_to_weather_event(
    injection_event::Vector{InjectionEvent},
    reservoir_properties::ReservoirProperties,
    domain::Domain3D
)::Vector{WeatherEvent}
    weather_events = [WeatherEvent(ie.timestamp, physical_volume_to_swim_volume(ie.injection_rate, reservoir_properties, domain)) for ie in injection_event]
    return weather_events
end