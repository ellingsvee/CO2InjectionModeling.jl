using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using Printf

# Setup similar to debugging_single_layer.jl
boundary_condition = :closed
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

# Load injection location from feeder chimney data
xy, (utm_x, utm_y, depth) = load_feeder_location(topography)
println("Injection location loaded from feeder data:")
println("  Grid indices: $xy")
println("  UTM coordinates: ($utm_x, $utm_y)")
println("  Depth: $depth m")
injection_events = generate_sleipner_injection_events(layers, xy)

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

# Generate multi-layer animations
println("\n=== Generating Multi-Layer Animations ===")
println("Generating height-based animation...")
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

println("\nGenerating saturation-based animation...")
animate_multi_layer_saturation(
    layers,
    seqs,
    leakage_states,
    domain,
    rprops;
    output_file="multi_layer_saturation.gif",
    num_frames=30,
    end_time=15.0,
    fps=3,
    colormap=:viridis
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

# Verify mass conservation
println("\n=== Mass Conservation Verification ===")
println("Checking that total_injected = total_stored + total_leaked at each snapshot...\n")

# Print table header
println("Time (yr) | Injected (m³) | Stored (m³) | Leaked (m³) | Error (m³) | Error (%)")
println("-" ^ 85)

# Track maximum error
max_error = 0.0
max_error_percent = 0.0

for snapshot in snapshots
    error_m3 = snapshot.mass_balance_error_m3
    error_pct = snapshot.mass_balance_error_percent

    global max_error = max(max_error, abs(error_m3))
    global max_error_percent = max(max_error_percent, error_pct)

    @printf("%8.2f  | %13.2e | %11.2e | %11.2e | %10.2e | %8.4f\n",
            snapshot.timestamp,
            snapshot.total_injected_m3,
            snapshot.total_stored_m3,
            snapshot.total_leaked_m3,
            error_m3,
            error_pct)
end

println("-" ^ 85)
println("\nMass Conservation Summary:")
println("  Maximum absolute error: $(max_error) m³")
println("  Maximum relative error: $(max_error_percent) %")

# Check if mass is conserved within tolerance
tolerance_percent = 0.01  # 0.01% tolerance
if max_error_percent < tolerance_percent
    println("  ✓ Mass is conserved (within $(tolerance_percent)% tolerance)")
else
    println("  ✗ Mass conservation violated! Error exceeds $(tolerance_percent)% tolerance")
end

println("\nTest completed!")
