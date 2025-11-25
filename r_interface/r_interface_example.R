# CO2 Injection Simulation Example
# Demonstrates the improved R interface for CO2InjectionModeling.jl

library(JuliaCall)

# Setup Julia and load package
julia_setup()

# REMEMBER TO SET WORKING DIR TO CO2InjectionModeling.jl
# Alternatively: Install the package globally if you want to call the lib from
# another project.
julia_command('using Pkg; Pkg.activate(".")')

# This line might take a bit of time...
julia_command('using CO2InjectionModeling')

# ============================================================================
# EXAMPLE 1: Simple simulation with Sleipner defaults
# ============================================================================

cat("=== Example 1: Simple Simulation with Sleipner Defaults ===\n\n")

# Step 1: Setup simulator
cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator", boundary_condition = "open") # Alternatively "closed" for closed BCs
print(setup_result)

nx <- setup_result$nx
ny <- setup_result$ny
n_layers <- setup_result$n_layers

# Step 2: Use Sleipner default reservoir properties
cat("\nUsing Sleipner default reservoir properties...\n")
config_result <- julia_call("setup_sleipner_reservoir")
print(config_result)

# Step 3: Setup injection scenario
cat("\nSetting up injection scenario...\n")

# Historical Sleipner injection rates (1996-2010)
rates_mt <- c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
              0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
              0.82, 0.86, 0.76)  # Mt/year

# Convert to m³/year (CO2 density at bottom layer = 570 kg/m³)
co2_density <- 570.0
rates_m3 <- rates_mt * 1e9 / co2_density

n_times <- length(rates_m3)

# Create injection matrix for layer 1 (bottom layer): n_times × nx × ny
layer1_injection <- array(0, dim = c(n_times, nx, ny))

# Inject at location (32, 59) - approximate center of Sleipner field
for (i in 1:n_times) {
  layer1_injection[i, 32, 59] <- rates_m3[i]
}

# Create zero injection for other layers (1 × nx × ny)
zero_injection <- array(0, dim = c(1, nx, ny))

# Build list of injection matrices (one per layer)
injection_matrices <- list(
  layer1_injection,  # Layer 1 (bottom)
  zero_injection,    # Layer 2
  zero_injection,    # Layer 3
  zero_injection,    # Layer 4
  zero_injection,    # Layer 5
  zero_injection,    # Layer 6
  zero_injection,    # Layer 7
  zero_injection,    # Layer 8
  zero_injection     # Layer 9 (top)
)

# Step 4: Run simulation
cat("\nRunning simulation (this may take a moment on first run)...\n")
sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 15.0,
                         time_step = 1.0,
                         injection_rate_matrices = injection_matrices,
                         verbose = FALSE)

# Step 5: Check results
if (sim_result$status == "success") {
  cat("\n=== Simulation Successful! ===\n")
  cat("Timepoints:", sim_result$timepoints, "\n")
  cat("Total CO2 volumes (m³):\n")
  print(sim_result$total_co2_volumes)
  cat("\nFinal total volume:", tail(sim_result$total_co2_volumes, 1), "m³\n")
} else {
  cat("\nSimulation failed:", sim_result$message, "\n")
  if (!is.null(sim_result$stacktrace)) {
    cat("Stacktrace:\n", sim_result$stacktrace, "\n")
  }
}


# ============================================================================
# EXAMPLE 2: Custom reservoir properties with multiple injection wells
# ============================================================================

cat("\n\n=== Example 2: Custom Properties with Multiple Wells ===\n\n")

# Setup simulator (reuse from above)
cat("Configuring custom reservoir properties...\n")

# Use layer-specific density differences
brine_density <- 1020
co2_densities <- c(570, 542.5, 515, 487.5, 460, 432.5, 405, 377.5, 350)
density_diffs <- brine_density - co2_densities

custom_config <- julia_call("configure_reservoir",
                             porosity = rep(0.35, 9),
                             residual_co2_sat = rep(0.25, 9),
                             irreducible_water_sat = rep(0.3, 9),
                             shale_pressure_threshold = rep(98000.0, 9),
                             brine_co2_density_diff = density_diffs,
                             residual_leakage_time = rep(1.0, 9),
                             layer_specific = TRUE)
print(custom_config)

# Create injection scenario with two wells
cat("\nSetting up two-well injection scenario...\n")

n_times <- 10
well1_rates <- seq(0.5, 1.0, length.out = n_times) * 1e9 / 570  # Ramping up
well2_rates <- rep(0.8, n_times) * 1e9 / 570  # Constant rate

# Layer 1 injection with two wells
layer1_injection_2wells <- array(0, dim = c(n_times, nx, ny))

for (i in 1:n_times) {
  layer1_injection_2wells[i, 32, 59] <- well1_rates[i]  # Well 1
  layer1_injection_2wells[i, 35, 62] <- well2_rates[i]  # Well 2
}

# Build injection matrices
injection_matrices_2wells <- list(
  layer1_injection_2wells,
  zero_injection, zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection, zero_injection, zero_injection
)

# Run simulation
cat("\nRunning two-well simulation...\n")
sim_result_2wells <- julia_call("run_simulation",
                                 start_time = 0.0,
                                 end_time = 10.0,
                                 time_step = 1.0,
                                 injection_rate_matrices = injection_matrices_2wells,
                                 verbose = FALSE)

if (sim_result_2wells$status == "success") {
  cat("\n=== Two-Well Simulation Successful! ===\n")
  cat("Final total volume:", tail(sim_result_2wells$total_co2_volumes, 1), "m³\n")
} else {
  cat("\nTwo-well simulation failed:", sim_result_2wells$message, "\n")
}

cat("\n=== All examples completed! ===\n")
