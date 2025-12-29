using CO2InjectionModeling
using SurfaceWaterIntegratedModeling

# Setup similar to debugging_single_layer.jl
boundary_condition = :closed
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

trap_topo = layers[1].trap_structure.topography
xy = CartesianIndex(div(size(trap_topo, 1), 2), div(size(trap_topo, 2), 2))
injection_events = generate_sleipner_injection_events(layers, xy)

# Create reservoir properties with a VERY LOW leakage height to force leakage
layer_idx = 1
rprops = generate_reservoir_properties_for_sleipner_layers()

# Run fill_layer with leakage detection
seqs, leakage_states = fill_layers(
    layers,
    domain,
    rprops,
    injection_events;
    verbose=false
);

# Generate multi-layer animation
println("\n=== Generating Multi-Layer Animation ===")
animate_multi_layer_filling(
    layers,
    seqs,
    domain;
    output_file="multi_layer_filling.gif",
    num_frames=30,
    end_time=15.0,
    fps=3,
    max_CO2_height=rprops[1].leakage_height
)

# Generate reservoir snapshots for analysis
println("\n=== Generating Reservoir Snapshots ===")
snapshots = generate_reservoir_snapshots(
    layers,
    seqs,
    leakage_states,
    domain,
    rprops,
    injection_events;
    num_snapshots=16,
    start_time=0.0,
    end_time=15.0,
    verbose=true
)

# Print summaries at key timepoints
println("\n=== Snapshot Summaries ===")
for snapshot in snapshots[1:4:end]  # Every 4th snapshot
    print_snapshot_summary(snapshot)
end

# Print detailed layer info for final snapshot
println("\n=== Final Snapshot - Layer Details ===")
final_snapshot = snapshots[end]
for layer_snapshot in final_snapshot.layer_snapshots
    print_layer_snapshot_summary(layer_snapshot)
end

# Generate timeseries plots
println("\n=== Generating Timeseries Plots ===")

# Plot CO2 volumes by layer (subplots)
println("Generating layer volumes plot...")
plot_layer_volumes_timeseries(
    snapshots;
    output_file="layer_volumes_timeseries.png",
    title="CO2 Volume by Layer Over Time"
)

# Plot CO2 fractions (stacked area chart)
println("Generating layer fractions plot...")
plot_layer_fractions_timeseries(
    snapshots;
    output_file="layer_fractions_timeseries.png",
    title="CO2 Distribution Across Layers Over Time"
)

println("\nTest completed!")
