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
xy = CartesianIndex(div(size(trap_topo, 1), 2), div(size(trap_topo, 2), 2))
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

# Generate animations
# Cross-section animation (all layers combined)
animate_trap_filling_multilayer(snapshots, layers, lithology, domain,
    output_file="multilayer_xsection.gif", direction=:x)

# Bird's eye grid animation (3x3 subplots)
animate_trap_filling_birdseye_multilayer(snapshots, layers, domain,
    output_file="multilayer_birdseye.gif")

# Plot total CO2 volume over time
plot_total_co2_volume(summary, output_file="co2_volume.png")

# Plot layer-wise CO2 volumes over time
plot_layer_co2_volumes(summary, output_file="layerwise_co2_volumes.png")