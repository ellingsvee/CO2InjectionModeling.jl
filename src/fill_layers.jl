using SurfaceWaterIntegratedModeling

export fill_layers

function fill_layers(
        layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::ReservoirProperties,
        injection_events::Vector{InjectionEvent};
        num_snapshots::Int=10,
        start_time::Float64=0.0,
        end_time::Float64=15.0,
        verbose::Bool=false

)
    weather_events = [WeatherEvent(ie.timestamp, physical_volume_to_swim_volume(ie.injection_rate, reservoir_properties, domain)) for ie in injection_events]

    n_layers = length(layers)
    seqs = Vector{Vector{SpillEvent}}(undef, n_layers)
    leakage_states = Vector{LeakageState}(undef, n_layers)

    weather_events_layer = weather_events
    for layer_idx in 1:n_layers
        if verbose
            println("Filling layer $layer_idx / $n_layers")
        end

        tstruct = layers[layer_idx].trap_structure

        seq, leakage_state = fill_layer_co2(
            tstruct,
            domain,
            reservoir_properties,
            weather_events_layer;
            verbose=verbose
        )

        seqs[layer_idx] = seq
        leakage_states[layer_idx] = leakage_state

        # TODO: Update the weather events for the next layer based on leakage from this layer
        # The leakage_state contains information about the leaked CO2 volume and timing
        
    end

    # TODO: For evey layer, create num_snapshots snapshots between start_time and end_time
    # Use the SimulationSnapshot and SimulationLayerSnapshot structs from types.jl as an outline, but possibly extend them with more fields as needed.
    # See the simulation_layer_snapshots_from_spill_events for some inspiration for how to do this. However, this is legacy code that needs to be adapted to work with multiple layers.
    # Recall that each layer has its own spill event sequence and leakage_events.


    return seqs, leakage_states

end

