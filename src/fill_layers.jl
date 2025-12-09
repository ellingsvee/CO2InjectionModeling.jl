using SurfaceWaterIntegratedModeling

export fill_layers, convert_injection_event_to_weather_event, create_next_layer_weather_events, generate_multi_level_snapshots





function fill_layers(
        layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}},
        injection_events::Vector{Vector{InjectionEvent}};
        num_snapshots::Int=10,
        start_time::Float64=0.0,
        end_time::Float64=15.0,
        verbose::Bool=false

)
    n_layers = length(layers)

    # If a single ReservoirProperties is provided, replicate it for all layers
    if reservoir_properties isa ReservoirProperties
        reservoir_properties = fill(reservoir_properties, n_layers)
    end


    seqs = Vector{Vector{SpillEvent}}(undef, n_layers)
    leakage_states = Vector{LeakageState}(undef, n_layers)

    weather_events_layer = convert_injection_event_to_weather_event(injection_events[1], reservoir_properties[1], domain)
    for layer_idx in 1:n_layers
        if verbose
            println("Filling layer $layer_idx / $n_layers")
        end

        tstruct = layers[layer_idx].trap_structure

        seq, leakage_state = fill_layer(
            tstruct,
            domain,
            reservoir_properties[layer_idx],
            weather_events_layer;
            verbose=verbose
        )

        seqs[layer_idx] = seq
        leakage_states[layer_idx] = leakage_state

        if layer_idx < n_layers
            weather_events_next_layer = convert_injection_event_to_weather_event(injection_events[layer_idx + 1], reservoir_properties[layer_idx + 1], domain)
            weather_events_layer = create_next_layer_weather_events(
                leakage_state,
                layers[layer_idx].trap_structure,
                layers[layer_idx + 1].trap_structure,
                weather_events_next_layer,
                weather_events_layer
            )

        end
    end

    # TODO: For evey layer, create num_snapshots snapshots between start_time and end_time
    # Use the SimulationSnapshot and SimulationLayerSnapshot structs from types.jl as an outline, but possibly extend them with more fields as needed.
    # See the simulation_layer_snapshots_from_spill_events for some inspiration for how to do this. However, this is legacy code that needs to be adapted to work with multiple layers.
    # Recall that each layer has its own spill event sequence and leakage_events.


    snapshots = generate_multi_level_snapshots(
        layers,
        seqs,
        domain,
        reservoir_properties,
        injection_events,
        num_snapshots,
        start_time,
        end_time
    )


    return seqs, leakage_states, snapshots

end

function convert_injection_event_to_weather_event(
        injection_event::Vector{InjectionEvent},
        reservoir_properties::ReservoirProperties,
        domain::Domain3D
    )::Vector{WeatherEvent}
    weather_events = [WeatherEvent(ie.timestamp, physical_volume_to_swim_volume(ie.injection_rate, reservoir_properties, domain)) for ie in injection_event]
    return weather_events
end

"""
Create weather events for the next layer based on leakage from the current layer.

When CO2 leaks from a trap in the current layer, it becomes injection into the next layer.
The leakage location (bottom of the leaked trap) becomes the injection location.

Simplified model (no residual trapping):
- All CO2 that leaks goes through immediately at the inflow rate
- Leakage starts at le.timestamp and continues at the inflow rate to that trap

Parameters:
- leakage_state: LeakageState from the current layer
- current_tstruct: TrapStructure of the current layer
- next_tstruct: TrapStructure of the next layer
- next_layer_base_weather: Base weather events for the next layer (from direct injection)
- current_layer_weather: Weather events used for the current layer (to get inflow rates)

Returns:
- Vector{WeatherEvent}: Combined weather events for the next layer
"""
function create_next_layer_weather_events(
    leakage_state::LeakageState,
    current_tstruct::TrapStructure,
    next_tstruct::TrapStructure,
    next_layer_base_weather::Vector{WeatherEvent},
    current_layer_weather::Vector{WeatherEvent}
)
    nx, ny = size(next_tstruct.topography)

    # If no leakage occurred, just return the base weather events for the next layer
    if isempty(leakage_state.leakage_events)
        if isempty(next_layer_base_weather)
            return [WeatherEvent(0.0, zeros(Float64, nx, ny))]
        end
        return next_layer_base_weather
    end

    # Collect all important timestamps:
    # 1. From base weather events for this layer
    # 2. From leakage events (when leakage starts)
    # 3. From current layer weather events (rate changes)
    all_timestamps = Set{Float64}()
    push!(all_timestamps, 0.0)

    for we in next_layer_base_weather
        push!(all_timestamps, we.timestamp)
    end

    for we in current_layer_weather
        push!(all_timestamps, we.timestamp)
    end

    for le in leakage_state.leakage_events
        push!(all_timestamps, le.timestamp)
    end

    sorted_timestamps = sort(collect(all_timestamps))
    weather_events = Vector{WeatherEvent}()

    for timestamp in sorted_timestamps
        # Start with base injection rate for this layer at this time
        rain_rate = zeros(Float64, nx, ny)

        # Add base injection from next_layer_base_weather
        for we in reverse(next_layer_base_weather)
            if we.timestamp <= timestamp
                if we.rain_rate isa Matrix
                    rain_rate .+= we.rain_rate
                end
                break
            end
        end

        # Add leakage contributions from the current layer to the next layer
        # This includes:
        # 1. Ongoing injection in the leaked regions that is now redirected through the leak
        #    (the injection that was stopped in the current layer must continue to the next layer)
        for le in leakage_state.leakage_events
            # Skip if this timestamp is before the leakage started
            if timestamp < le.timestamp
                continue
            end

            # Find the leakage location in the current layer (where CO2 exits to next layer)
            leakage_location = find_leakage_location(le.trap_id, current_tstruct)
            injection_i, injection_j = leakage_location.I

            # Find the weather event active at the current timestamp
            # NOTE: current_layer_weather contains the ORIGINAL injection rates (not modified by leakage)
            # because fill_layer makes copies before applying leakage
            current_layer_weather_event = nothing
            for we in reverse(current_layer_weather)
                if we.timestamp <= timestamp
                    current_layer_weather_event = we
                    break
                end
            end
            @assert !isnothing(current_layer_weather_event) "No weather event found for current layer at time $timestamp"

            # Calculate the injection rate that should be redirected to the next layer
            # This is the injection that was happening in the source regions at this timestamp
            redirected_injection_rate = 0.0
            if !isempty(le.source_regions)
                for source_region in le.source_regions
                    region_mask = current_tstruct.regions .== source_region
                    # Get the injection rate in this source region at the current time
                    # This injection was stopped in the current layer (due to leakage) and
                    # must be redirected to the next layer
                    redirected_injection_rate += sum(current_layer_weather_event.rain_rate[region_mask])
                end
            end

            # Add the redirected injection to the leakage point in the next layer
            rain_rate[injection_i, injection_j] += redirected_injection_rate
        end

        push!(weather_events, WeatherEvent(timestamp, rain_rate))
    end

    return weather_events
end

function generate_multi_level_snapshots(
        layers::Vector{Layer},
        seqs::Vector{Vector{SpillEvent}},
        domain::Domain3D,
        reservoir_properties::Vector{ReservoirProperties},
        injection_events::Vector{Vector{InjectionEvent}},
        num_snapshots::Int,
        start_time::Float64,
        end_time::Float64
    )
    n_layers = length(layers);
    layer_snapshots = Vector{Vector{SimulationLayerSnapshot}}(undef, n_layers);


    for layer_idx in 1:n_layers
        tstruct = layers[layer_idx].trap_structure;
        seq = seqs[layer_idx];

        layer_snapshots[layer_idx] = simulation_layer_snapshots_from_spill_events(
            layers[layer_idx],
            seq,
            domain,
            reservoir_properties[layer_idx],
            injection_events[layer_idx];
            num_snapshots=num_snapshots,
            start_time=start_time,
            end_time=end_time
        );
    end

    # Combine layer snapshots into multi-level snapshots
    multi_level_snapshots = Vector{SimulationSnapshot}(undef, num_snapshots)
    for snapshot_idx in 1:num_snapshots
        layer_snaps = Vector{SimulationLayerSnapshot}(undef, n_layers)
        for layer_idx in 1:n_layers
            layer_snaps[layer_idx] = layer_snapshots[layer_idx][snapshot_idx]
        end

        multi_level_snapshots[snapshot_idx] = SimulationSnapshot(
            layer_snaps[1].timestamp,
            sum(ls.injected_volume for ls in layer_snaps),
            sum(ls.co2_volume for ls in layer_snaps),
            sum(ls.mobile_co2_volume for ls in layer_snaps),
            sum(ls.residual_trapped_co2_volume for ls in layer_snaps),
            layer_snaps
        )
    end

    return multi_level_snapshots
end