using CO2InjectionModeling
using SurfaceWaterIntegratedModeling

# Setup (same as run.jl)
boundary_condition = :closed
topography = load_sleipner_topography();
domain = create_domain_from_topography(topography, 1.0);
lithology = reconstruct_3d_lithology(topography, domain);
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition);

# Create CO2 saturation field and fill a trap
layer = 1;

# Find the trap that has the largest footprint
max_size = 0
max_trap = 0
for (trap_id, footprint) in enumerate(layers[layer].trap_structure.footprints)
    global max_size, max_trap
    if length(footprint) > max_size
        max_size = length(footprint)
        max_trap = trap_id
    end
end
trap = max_trap

# Get the topography dimensions from the trap structure (not domain!)
trap_topo = layers[layer].trap_structure.topography
xy = CartesianIndices(size(trap_topo))[layers[layer].trap_structure.footprints[trap][40]]

reservoir_props = ReservoirProperties()

injection_rate = fill(0.0, size(layers[layer].trap_structure.topography))
injection_rate[xy] = 1e6

injection_rate_2 = fill(0.0, size(layers[layer].trap_structure.topography))
injection_rate_2[xy[1]-25, xy[2] + 70] = 0.5e6

injection_event_bottom_layer = [InjectionEvent(0.0, injection_rate)]
injection_event_remaining_layers = [InjectionEvent(0.0, zeros(size(layers[layer].trap_structure.topography)))]

injection_events = Vector{Vector{InjectionEvent}}(undef, length(layers))
injection_events[1] = injection_event_bottom_layer
for i in 2:length(layers)
    if i == 4
        injection_events[i] = [InjectionEvent(0.0, injection_rate_2)]
    else
        injection_events[i] = injection_event_remaining_layers
    end
end

# Run simulation
seqs, leakage_states, snapshots = fill_layers(
    layers,
    domain,
    reservoir_props,
    injection_events;
    num_snapshots=14,
    start_time=0.0,
    end_time=15.0,
    verbose=false
)

println("Simulation completed with ", length(snapshots), " snapshots")

# Generate summary
println("\nGenerating simulation summary...")
summary = generate_simulation_summary(snapshots, layers, seqs)

println("\n=== Simulation Summary ===")
println("Number of layers: ", summary.num_layers)
println("Number of traps per layer: ", summary.num_traps_per_layer)
println("\nTimepoints: ", summary.timepoints)
println("\nTotal CO2 volumes: ", summary.total_co2_volumes)

println("\n--- Layer 1 Analysis ---")
println("Layer 1 volumes over time: ", summary.layer_co2_volumes[:, 1])

# Find the trap with most CO2 at the end
final_trap_volumes = summary.trap_co2_volumes[1][end, :]
max_trap_vol, max_trap_idx = findmax(final_trap_volumes)
println("\nTrap with most CO2 in Layer 1 at final time: Trap ", max_trap_idx)
println("  Volume: ", max_trap_vol)
println("  Percentage of total: ", summary.trap_co2_percentages[1][end, max_trap_idx], "%")

println("\n--- Layer 4 Analysis ---")
if summary.num_layers >= 4
    println("Layer 4 volumes over time: ", summary.layer_co2_volumes[:, 4])
    final_trap_volumes_4 = summary.trap_co2_volumes[4][end, :]
    if any(final_trap_volumes_4 .> 0)
        max_trap_vol_4, max_trap_idx_4 = findmax(final_trap_volumes_4)
        println("\nTrap with most CO2 in Layer 4 at final time: Trap ", max_trap_idx_4)
        println("  Volume: ", max_trap_vol_4)
        println("  Percentage of total: ", summary.trap_co2_percentages[4][end, max_trap_idx_4], "%")
    else
        println("No CO2 in Layer 4")
    end
end

println("\n✓ Summary generation successful!")
