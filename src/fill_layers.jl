using SurfaceWaterIntegratedModeling
export fill_layers

"""
    fill_layers(layers, domain, reservoir_properties, injection_events; verbose=false)
        -> (seqs, leakage_states, weather_events_per_layer)

Simulate CO2 migration through a stack of geological layers, propagating
leakage upward from each layer to the next.

Layers are processed in index order: `layers[1]` (deepest) is filled first
using `injection_events[1]`.  CO2 that leaks through the caprock of layer `i`
is converted to [`InjectionEvent`](@ref)-equivalent `WeatherEvent`s and drives
filling in layer `i+1`, and so on.

# Arguments
- `layers`: `Vector{Layer}` from [`analyze_base_surfaces`](@ref), deepest first
- `domain`: [`Domain3D`](@ref) from [`create_domain`](@ref)
- `reservoir_properties`: A single [`ReservoirProperties`](@ref) (applied to
  all layers) or a `Vector{ReservoirProperties}` (one per layer)
- `injection_events`: `Vector{Vector{InjectionEvent}}` — one inner vector per
  layer.  Layers without direct injection receive a single zero-rate event.
- `verbose`: Print per-layer progress (default `false`)

# Returns
- `seqs`: `Vector{Vector{SpillEvent}}` — fill sequence per layer
- `leakage_states`: `Vector{LeakageState}` — leakage tracking per layer
- `weather_events_per_layer`: `Vector{Vector{WeatherEvent}}` — effective
  injection (including leakage from below) per layer

# Example
```julia
seqs, leakage_states, weather_events_per_layer = fill_layers(
    layers, domain, [rp, rp, rp_caprock], injection_events)
```
"""
function fill_layers(
    layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties,Vector{ReservoirProperties}},
    injection_events::Vector{Vector{InjectionEvent}};
    verbose::Bool=false,
    leakage_radius::Int=0,
)
    n_layers = length(layers)

    seqs = Vector{Vector{SpillEvent}}(undef, n_layers)
    leakage_states = Vector{LeakageState}(undef, n_layers)
    weather_events_per_layer = Vector{Vector{WeatherEvent}}(undef, n_layers)

    # In the case where each layer has the same properties
    if reservoir_properties isa ReservoirProperties
        reservoir_properties = fill(reservoir_properties, n_layers)
    end

    weather_events_layer = convert_injection_event_to_weather_event(injection_events[1], reservoir_properties[1], domain)

    for layer_idx in 1:n_layers
        verbose && println("Filling layer $layer_idx / $n_layers")

        weather_events_per_layer[layer_idx] = weather_events_layer

        tstruct = layers[layer_idx].trap_structure

        seq, leakage_state = fill_sequence_with_leakage(
            tstruct,
            reservoir_properties[layer_idx],
            weather_events_layer;
            verbose=verbose
        )

        seqs[layer_idx] = seq
        leakage_states[layer_idx] = leakage_state


        if layer_idx < n_layers
            weather_events_next_layer = convert_injection_event_to_weather_event(injection_events[layer_idx+1], reservoir_properties[layer_idx+1], domain)
            weather_events_layer = generate_leakage_weather_events(
                seq,
                leakage_state,
                tstruct,
                reservoir_properties[layer_idx],
                reservoir_properties[layer_idx+1],
                weather_events_next_layer;
                leakage_radius=leakage_radius,
                target_regions=layers[layer_idx+1].trap_structure.regions,
            )

        end
    end

    return seqs, leakage_states, weather_events_per_layer
end