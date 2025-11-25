using CO2InjectionModeling
using SurfaceWaterIntegratedModeling

boundary_condition = :closed
topography = load_sleipner_topography();
domain = create_domain_from_topography(topography, 1.0);
lithology = reconstruct_3d_lithology(topography, domain);
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition);


# Generate reservoir properties
reservoir_properties = generate_reservoir_properties_for_sleipner_layers()

# Define injection events
# Make xy in the center of the domain
trap_topo = layers[1].trap_structure.topography
xy = CartesianIndex(div(size(trap_topo, 1), 2), div(size(trap_topo, 2), 4))
injection_events = generate_sleipner_injection_events(layers, xy)



# Run injection in a SINGLE layer for debugging
layer_idx = 1
weather_events_layer = convert_injection_event_to_weather_event(injection_events[layer_idx], reservoir_properties[layer_idx], domain)
tstruct = layers[layer_idx].trap_structure


seq, leakage_state = fill_layer(
    tstruct,
    domain,
    reservoir_properties[layer_idx],
    # [weather_events_layer[1]]; # Using just a single weater event (not varying injection rate over time works)
    weather_events_layer;
    verbose=false
)

weather_events_next_layer = convert_injection_event_to_weather_event(injection_events[layer_idx + 1], reservoir_properties[layer_idx + 1], domain)
weather_events_layer_updated = create_next_layer_weather_events(
    leakage_state,
    layers[layer_idx].trap_structure,
    layers[layer_idx + 1].trap_structure,
    weather_events_next_layer,
    weather_events_layer
)

# Inspect results
println(leakage_state.first_leakage_time) # 1.8697300698015689
println(length(leakage_state.leakage_events)) # 5

# Amount at first leakage time
tstates = trap_states_at_timepoints(tstruct, seq, [leakage_state.first_leakage_time])
amounts_at_leakage = [e[2] for e in tstates]
total_amounts_at_leakage = sum(amounts_at_leakage[1])

# Amount 5 time units after first leakage
tstates = trap_states_at_timepoints(tstruct, seq, [leakage_state.first_leakage_time + 5])
amounts_time_after_leakage_leakage = [e[2] for e in tstates]
total_amounts_time_after_leakage_leakage = sum(amounts_time_after_leakage_leakage[1])

println("Total amount at leakage time: ", total_amounts_at_leakage)
println("Total amount 5 time units after leakage: ", total_amounts_time_after_leakage_leakage)