# Simple CO2 Injection Simulation Example
# Minimal example to get started quickly

library(JuliaCall)

# Setup Julia and load package
julia_setup()
julia_command('using Pkg; Pkg.activate(".")')
julia_command('using CO2InjectionModeling')
julia_command('using .CO2RInterface')

cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator",
                           boundary_condition = "open")
print(setup_result)

cat("\nConfiguring reservoir properties...\n")
config_result <- julia_call("configure_reservoir",
                            porosity = 0.4,
                            residual_co2_sat = 0.2,
                            irreducible_water_sat = 0.3,
                            shale_pressure_threshold = 98000.0,
                            brine_co2_density_diff = 450.0,
                            residual_leakage_time = 1.0,
                            layer_specific = FALSE)
print(config_result)

cat("\nRunning simulation...\n")

# Simple injection: constant rate at one location for 5 years
injection_times <- c(0, 1, 2, 3, 4)
injection_rate <- 1e6  # m³/year
injection_i <- rep(32L, 5)
injection_j <- rep(59L, 5)
injection_amounts <- rep(injection_rate, 5)
injection_layers <- rep(1L, 5)  # Bottom layer

sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 5.0,
                         time_step = 1.0,
                         injection_times = injection_times,
                         injection_locations_i = injection_i,
                         injection_locations_j = injection_j,
                         injection_amounts = injection_amounts,
                         injection_layer_indices = injection_layers,
                         num_snapshots = 5L,
                         verbose = TRUE)

if (sim_result$status == "success") {
  cat("\nSimulation successful!\n")
  cat("Timepoints:", sim_result$timepoints, "\n")
  cat("Total volumes:", sim_result$total_co2_volumes, "\n")
} else {
  cat("\nSimulation failed:", sim_result$message, "\n")
}
