# Unit Conversion

CO2BatchFill uses SWIM's internal dimensionless volume units internally.
These functions convert between SWIM units and physical volumes and handle
the conversion of injection events to SWIM weather events.

## Volume conversion

```@docs
CO2BatchFill.swim_volume_to_physical_volume
CO2BatchFill.physical_volume_to_swim_volume
CO2BatchFill.full_volume_to_rock_volume_scaling
CO2BatchFill.unit_volume_to_physical_scaling
```

## Injection event conversion

```@docs
CO2BatchFill.convert_injection_event_to_weather_event
```
