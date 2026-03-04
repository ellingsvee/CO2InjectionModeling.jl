using SurfaceWaterIntegratedModeling
export fill_layers

"""
Extend the filling to multiple layers, accounting for leakage between layers.
"""
function fill_layers(
    layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties,Vector{ReservoirProperties}},
    injection_events::Vector{Vector{InjectionEvent}};
    verbose::Bool=false
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
                weather_events_next_layer
            )

        end
    end

    return seqs, leakage_states, weather_events_per_layer
end