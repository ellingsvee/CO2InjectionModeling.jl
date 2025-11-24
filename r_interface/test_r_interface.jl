# Test script to verify R interface functions work correctly
# Run this from Julia to test the interface without needing R

using CO2InjectionModeling

println("=" ^ 60)
println("Testing CO2RInterface Functions")
println("=" ^ 60)

# Test 1: setup_simulator
println("\n[Test 1] setup_simulator")
println("-" ^ 40)
setup_result = setup_simulator(
    data_path="sleipner/depth_surfaces/",
    boundary_condition="open"
)
println("Status: ", setup_result["status"])
if setup_result["status"] == "success"
    println("  n_layers: ", setup_result["n_layers"])
    println("  nx: ", setup_result["nx"])
    println("  ny: ", setup_result["ny"])
    println("  boundary_condition: ", setup_result["boundary_condition"])
    println("✓ Test passed")
else
    println("✗ Test failed: ", setup_result["message"])
    exit(1)
end

# Test 2: configure_reservoir (uniform properties)
println("\n[Test 2] configure_reservoir (uniform)")
println("-" ^ 40)
config_result = configure_reservoir(
    0.4,    # porosity
    0.2,    # residual_co2_sat
    0.3,    # irreducible_water_sat
    98000.0, # shale_pressure_threshold
    450.0,  # brine_co2_density_diff
    1.0;    # residual_leakage_time
    layer_specific=false
)
println("Status: ", config_result["status"])
if config_result["status"] == "success"
    println("  n_layers: ", config_result["n_layers"])
    println("✓ Test passed")
else
    println("✗ Test failed: ", config_result["message"])
    exit(1)
end

# Test 3: configure_reservoir (layer-specific properties)
println("\n[Test 3] configure_reservoir (layer-specific)")
println("-" ^ 40)
n_layers = setup_result["n_layers"]
config_result2 = configure_reservoir(
    fill(0.4, n_layers),
    fill(0.2, n_layers),
    fill(0.3, n_layers),
    fill(98000.0, n_layers),
    [450.0, 477.5, 505.0, 532.5, 560.0, 587.5, 615.0, 642.5, 670.0],
    fill(1.0, n_layers);
    layer_specific=true
)
println("Status: ", config_result2["status"])
if config_result2["status"] == "success"
    println("  n_layers: ", config_result2["n_layers"])
    println("✓ Test passed")
else
    println("✗ Test failed: ", config_result2["message"])
    exit(1)
end

# Test 4: run_simulation (small test)
println("\n[Test 4] run_simulation (small)")
println("-" ^ 40)

# Simple injection: constant rate for 3 years
injection_times = [0.0, 1.0, 2.0]
injection_i = [32, 32, 32]
injection_j = [59, 59, 59]
injection_amounts = [1e6, 1e6, 1e6]
injection_layers = [1, 1, 1]

sim_result = run_simulation(
    0.0,    # start_time
    3.0,    # end_time
    1.0,    # time_step
    injection_times,
    injection_i,
    injection_j,
    injection_amounts,
    injection_layers;
    num_snapshots=3,
    verbose=true
)

println("\nStatus: ", sim_result["status"])
if sim_result["status"] == "success"
    println("  Number of snapshots: ", length(sim_result["timepoints"]))
    println("  Timepoints: ", sim_result["timepoints"])
    println("  Total volumes: ", sim_result["total_co2_volumes"])
    println("  Number of layers: ", sim_result["num_layers"])
    println("  Number of traps per layer: ", sim_result["num_traps_per_layer"])

    # Verify data dimensions
    @assert length(sim_result["timepoints"]) == 3
    @assert length(sim_result["total_co2_volumes"]) == 3
    @assert size(sim_result["layer_co2_volumes"], 1) == 3
    @assert size(sim_result["layer_co2_volumes"], 2) == sim_result["num_layers"]
    @assert length(sim_result["trap_co2_volumes"]) == sim_result["num_layers"]
    @assert length(sim_result["trap_co2_percentages"]) == sim_result["num_layers"]

    println("✓ Test passed")
else
    println("✗ Test failed: ", sim_result["message"])
    if haskey(sim_result, "stacktrace")
        println("\nStacktrace:")
        println(sim_result["stacktrace"])
    end
    exit(1)
end

# Test 5: Error handling - call run_simulation with missing parameters
println("\n[Test 5] Error handling")
println("-" ^ 40)

# Test with mismatched vector lengths
error_result = run_simulation(
    0.0, 3.0, 1.0,
    [0.0, 1.0],  # 2 times
    [32],  # but only 1 location
    [59],
    [1e6],
    [1];
    num_snapshots=1, verbose=false
)

println("Status: ", error_result["status"])
if error_result["status"] == "error"
    println("  Message: ", error_result["message"])
    println("✓ Test passed (error caught correctly)")
else
    println("✗ Test failed: Expected error but got success")
    exit(1)
end

# Test 6: Full workflow test (like Sleipner historical)
println("\n[Test 6] Full workflow (Sleipner-like)")
println("-" ^ 40)

# Re-setup
setup_result = setup_simulator(boundary_condition="closed")
println("Setup: ", setup_result["status"])

# Configure with default Sleipner properties
config_result = configure_reservoir(
    0.4, 0.2, 0.3, 98000.0, 450.0, 1.0;
    layer_specific=false
)
println("Configure: ", config_result["status"])

# Run with simplified Sleipner injection pattern
annual_rates_mt = [0.07, 0.67, 0.85, 0.94, 0.94]
co2_density = 570.0
injection_rates_m3 = annual_rates_mt .* 1e9 ./ co2_density

sim_result = run_simulation(
    0.0,
    5.0,
    1.0,
    Float64.(0:4),
    fill(32, 5),
    fill(59, 5),
    injection_rates_m3,
    fill(1, 5);
    num_snapshots=5,
    verbose=false
)

println("Simulation: ", sim_result["status"])
if sim_result["status"] == "success"
    final_volume = sim_result["total_co2_volumes"][end]
    total_injected = sum(injection_rates_m3)
    println("  Total injected: $(total_injected) m³")
    println("  Final volume: $(final_volume) m³")
    println("  Ratio: $(final_volume / total_injected)")

    # Check that we stored some CO2
    @assert final_volume > 0 "No CO2 was stored"
    @assert sim_result["num_layers"] == 9 "Expected 9 layers for Sleipner"

    println("✓ Test passed")
else
    println("✗ Test failed: ", sim_result["message"])
    exit(1)
end

println("\n" * "=" ^ 60)
println("All tests passed! ✓")
println("=" ^ 60)
println("\nThe R interface is ready to use.")
println("See examples/R_INTERFACE_README.md for usage instructions.")
