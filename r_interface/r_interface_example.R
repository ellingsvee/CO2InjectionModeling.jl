# Simple CO2 Injection Simulation Example
# Minimal example to get started quickly

library(JuliaCall)

# Setup Julia and load package
julia_setup()

# REMEMBER TO SET WORKING DIR TO CO2InjectionModeling.jl
# Alternatively: Install the package globally if you want to call the lib from
# another project.
julia_command('using Pkg; Pkg.activate(".")')

# These lines might take a bit of time...
julia_command('using CO2InjectionModeling')

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

# Define injection scenario: 15 years in bottom layer (L1)
# Years 0-4: 0.5 Mt/year
# Years 5-14: 1.0 Mt/year
n_years <- 15
injection_times <- as.numeric(seq(0, n_years - 1))  # Must be Float64

# Convert from Mt/year to m³/year (CO2 density = 570 kg/m³)
co2_density <- 570.0  # kg/m³
rates_mt <- c(rep(0.5, n_years))  # Mt/year
# rates_mt <- c(rep(0.5, 5), rep(1.0, 10)) # How to change the rate during the injection period.
injection_amounts <- as.numeric(rates_mt * 1e9 / co2_density)  # Must be Float64

# Injection location (grid coordinates, must be integers)
injection_i <- as.integer(rep(32, n_years))
injection_j <- as.integer(rep(59, n_years))

# All injection in bottom layer (layer 1, must be integers)
injection_layers <- as.integer(rep(1, n_years))

# The first time you run the simulation, it might take a little time.
# Due to the JIT in Julia, it will afterwards run very fast.
sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 15.0,
                         time_step = 1.0,
                         injection_times = injection_times,
                         injection_locations_i = injection_i,
                         injection_locations_j = injection_j,
                         injection_amounts = injection_amounts,
                         injection_layer_indices = injection_layers,
                         num_snapshots = 15L,
                         verbose = FALSE)

if (sim_result$status == "success") {
  cat("\nSimulation successful!\n")
  cat("Timepoints:", sim_result$timepoints, "\n")
  cat("Total volumes:", sim_result$total_co2_volumes, "\n")
} else {
  cat("\nSimulation failed:", sim_result$message, "\n")
}
