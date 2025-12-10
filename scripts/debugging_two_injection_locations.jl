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
include("debugging_utils.jl")
injection_events = generate_injection_events(layers)


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


idxs = 1:length(snapshots)
for idx in idxs
    println("Snapshot $(idx): Time = $(snapshots[idx].timestamp) years")
    println("  Total injected volume: $(snapshots[idx].total_injected_volume) m³")
    println("  Total CO2 volume: $(snapshots[idx].total_co2_volume) m³")
end

# Bird's eye grid animation (3x3 subplots)
animate_trap_filling_birdseye_multilayer(snapshots, layers, domain,
    output_file="multilayer_birdseye.gif")
