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
trap_topo = layers[layer].trap_structure.topography
xy = CartesianIndex(div(size(trap_topo, 1), 2), div(size(trap_topo, 2), 4))

# injection_rate = fill(0.0, size(layers[layer].trap_structure.topography))
# injection_rate[xy] = 1e6
# injection_event_bottom_layer = [InjectionEvent(0.0, injection_rate)]
# injection_event_remaining_layers = [InjectionEvent(0.0, zeros(size(layers[layer].trap_structure.topography)))]

# injection_events = Vector{Vector{InjectionEvent}}(undef, length(layers))
# injection_events[1] = injection_event_bottom_layer
# for i in 2:length(layers)
#     injection_events[i] = injection_event_remaining_layers
# end

injection_events = generate_sleipner_injection_events(layers, xy)

# Run simulation
seqs, leakage_states, snapshots = fill_layers(
    layers,
    domain,
    reservoir_properties,
    injection_events;
    num_snapshots=14,
    start_time=0.0,
    end_time=15.0,
    verbose=false
);

println("Simulation completed with ", length(snapshots), " snapshots")

# Generate summary
println("\nGenerating simulation summary...")
summary = generate_simulation_summary(snapshots, layers, seqs);

println("\n=== Simulation Summary ===")
println("Number of layers: ", summary.num_layers)
println("Number of traps per layer: ", summary.num_traps_per_layer)
println("\nTimepoints: ", summary.timepoints)
println("\nTotal CO2 volumes: ", summary.total_co2_volumes)

# Generate animations
# Cross-section animation (all layers combined)
animate_trap_filling_multilayer(snapshots, layers, lithology, domain,
    output_file="multilayer_xsection.gif", direction=:x)

# Bird's eye grid animation (3x3 subplots)
animate_trap_filling_birdseye_multilayer(snapshots, layers, domain,
    output_file="multilayer_birdseye.gif")

