using CO2InjectionModeling
using SurfaceWaterIntegratedModeling

boundary_condition = :closed
topography = load_sleipner_topography();
domain = create_domain_from_topography(topography, 1.0);
lithology = reconstruct_3d_lithology(topography, domain);
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition);


# Generate reservoir properties
reservoir_properties = generate_reservoir_properties_for_sleipner_layers()



include("debugging_utils.jl")
injection_events = generate_injection_events(layers)





layer_idx = 1
tstruct = layers[layer_idx].trap_structure

n_layers = length(layers)
weather_events_layer = convert_injection_event_to_weather_event(injection_events[layer_idx], reservoir_properties[layer_idx], domain)


seq, leakage_state = fill_layer(
    tstruct,
    domain,
    reservoir_properties[layer_idx],
    weather_events_layer;
)


layer_snapshot = CO2InjectionModeling.simulation_layer_snapshots_from_spill_events(
    layers[layer_idx],
    seq,
    domain,
    reservoir_properties[layer_idx],
    injection_events[layer_idx];
    num_snapshots=14,
    start_time=0.0,
    end_time=15.0,  
);



idxs = 1:length(layer_snapshot)
for idx in idxs
    println("Snapshot $(idx): Time = $(layer_snapshot[idx].timestamp) years")
    println("  Total injected volume: $(layer_snapshot[idx].injected_volume) m³")
    println("  Total CO2 volume: $(layer_snapshot[idx].co2_volume) m³")
end

