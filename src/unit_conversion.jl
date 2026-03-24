using SurfaceWaterIntegratedModeling

export swim_volume_to_physical_volume, physical_volume_to_swim_volume
export full_volume_to_rock_volume_scaling, unit_volume_to_physical_scaling
export convert_injection_event_to_weather_event
export swim_volume_to_mass, physical_volume_to_mass

"""
    unit_volume_to_physical_scaling(domain) -> Float64

Return the factor that converts a SWIM "unit volume" (one vertical unit per
horizontal cell) to physical volume in m^3.  Equals `dx * dy`.
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

Convert SWIM internal volume units to physical volume in m^3.

SWIM volumes are dimensionless (height × cell count); this function applies
the pore-volume scaling (`porosity × (1 - Swi)`) and the horizontal cell area
(`dx × dy`) to obtain a physically meaningful volume.

# Arguments
- `volume`: Scalar or array in SWIM units
- `reservoir_properties`: [`ReservoirProperties`](@ref)
- `domain`: [`Domain3D`](@ref)

# Example
```julia
scale = swim_volume_to_physical_volume(1.0, rp, domain)  # m^3 per SWIM unit
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

Convert physical volume in m^3 to SWIM internal volume units.
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
    swim_volume_to_mass(volume, reservoir_properties, domain) -> Float64

Convert SWIM internal volume units to mass in kg.

Applies pore-volume scaling, cell area, and CO2 density to convert from
SWIM's dimensionless volume to physical mass.

# Arguments
- `volume`: Scalar or array in SWIM units
- `reservoir_properties`: [`ReservoirProperties`](@ref)
- `domain`: [`Domain3D`](@ref)
"""
function swim_volume_to_mass(
    volume::Union{AbstractArray{<:Real},Real},
    reservoir_properties::ReservoirProperties,
    domain::Domain3D
)
    physical_vol = swim_volume_to_physical_volume(volume, reservoir_properties, domain)
    return physical_volume_to_mass(physical_vol, reservoir_properties.co2_density)
end

"""
    physical_volume_to_mass(volume_m3, co2_density) -> Float64

Convert physical volume in m^3 to mass in kg.

# Arguments
- `volume_m3`: Volume in m^3 (scalar or array)
- `co2_density`: CO2 density in kg/m^3
"""
function physical_volume_to_mass(
    volume_m3::Union{AbstractArray{<:Real},Real},
    co2_density::Real
)
    return volume_m3 .* co2_density
end

"""
    convert_injection_event_to_weather_event(injection_events, reservoir_properties, domain)
        -> Vector{WeatherEvent}

Convert a vector of [`InjectionEvent`](@ref)s (physical m^3/yr rates) into SWIM
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