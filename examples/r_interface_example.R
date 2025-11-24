# CO2 Injection Modeling - R Interface Example
# This script demonstrates how to use the Julia CO2 injection simulator from R

# Install JuliaCall if needed
# install.packages("JuliaCall")

library(JuliaCall)

# Initialize Julia
# You may need to specify the Julia path explicitly, e.g.:
# julia_setup(JULIA_HOME = "/path/to/julia/bin")
julia_setup()

# Load the CO2InjectionModeling package
# Make sure Julia can find the package - you may need to:
# 1. Add the package to your Julia environment, or
# 2. Use julia_command to activate the project environment
julia_command('using Pkg; Pkg.activate(".")')
julia_command('using CO2InjectionModeling')

# Access the interface functions from the CO2RInterface module
julia_command('using .CO2RInterface')

cat("=== CO2 Injection Simulation Example ===\n\n")

# ===========================
# Step 1: Setup the simulator
# ===========================
cat("Step 1: Setting up simulator with Sleipner topography...\n")

setup_result <- julia_call("setup_simulator",
                           data_path = "sleipner/depth_surfaces/",
                           boundary_condition = "open")

if (setup_result$status == "success") {
  cat(sprintf("  Success! Loaded %d layers with grid size %d x %d\n",
              setup_result$n_layers,
              setup_result$nx,
              setup_result$ny))
  cat(sprintf("  Boundary condition: %s\n\n", setup_result$boundary_condition))
} else {
  stop(sprintf("Setup failed: %s", setup_result$message))
}

# ===========================
# Step 2: Configure reservoir properties
# ===========================
cat("Step 2: Configuring reservoir properties...\n")

# Option A: Use the same properties for all layers (simple)
config_result <- julia_call("configure_reservoir",
                            porosity = 0.4,
                            residual_co2_sat = 0.2,
                            irreducible_water_sat = 0.3,
                            shale_pressure_threshold = 98000.0,
                            brine_co2_density_diff = 450.0,
                            residual_leakage_time = 1.0,
                            layer_specific = FALSE)

# Option B: Use layer-specific properties (more realistic)
# The Sleipner field has 9 layers (L1-L9)
# CO2 density decreases with depth, so density difference increases
# co2_densities <- c(570.0, 542.5, 515.0, 487.5, 460.0, 432.5, 405.0, 377.5, 350.0)
# brine_density <- 1020.0
# density_diffs <- brine_density - co2_densities
#
# config_result <- julia_call("configure_reservoir",
#                             porosity = rep(0.4, 9),
#                             residual_co2_sat = rep(0.2, 9),
#                             irreducible_water_sat = rep(0.3, 9),
#                             shale_pressure_threshold = rep(98000.0, 9),
#                             brine_co2_density_diff = density_diffs,
#                             residual_leakage_time = rep(1.0, 9),
#                             layer_specific = TRUE)

if (config_result$status == "success") {
  cat(sprintf("  Success! Configured properties for %d layers\n\n",
              config_result$n_layers))
} else {
  stop(sprintf("Configuration failed: %s", config_result$message))
}

# ===========================
# Step 3: Run the simulation
# ===========================
cat("Step 3: Running simulation...\n")

# Define injection parameters based on Sleipner historical data (1996-2010)
# Annual injection rates in Mt/year
annual_rates_mt <- c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
                     0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
                     0.82, 0.86, 0.76)

# Convert Mt/year to m³/year using CO2 density at bottom layer (570 kg/m³)
co2_density_l1 <- 570.0  # kg/m³
injection_rates_m3 <- annual_rates_mt * 1e9 / co2_density_l1

# Injection times (years 0-14)
injection_times <- seq(0, 14)

# Injection location (approximately center of domain)
# These are 1-based indices for Julia
injection_i <- rep(32L, 15)  # i-index
injection_j <- rep(59L, 15)  # j-index

# Inject into bottom layer (layer 1)
injection_layers <- rep(1L, 15)

# Run the simulation
cat("  Running simulation from year 0 to 15...\n")
sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 15.0,
                         time_step = 1.0,
                         injection_times = injection_times,
                         injection_locations_i = injection_i,
                         injection_locations_j = injection_j,
                         injection_amounts = injection_rates_m3,
                         injection_layer_indices = injection_layers,
                         num_snapshots = 14L,
                         verbose = TRUE)

if (sim_result$status != "success") {
  stop(sprintf("Simulation failed: %s\n%s",
               sim_result$message,
               sim_result$stacktrace))
}

cat("\n  Simulation completed successfully!\n\n")

# ===========================
# Step 4: Analyze results
# ===========================
cat("=== Simulation Results ===\n\n")

# Extract results
timepoints <- sim_result$timepoints
total_volumes <- sim_result$total_co2_volumes
layer_volumes <- sim_result$layer_co2_volumes
n_layers <- sim_result$num_layers
n_traps_per_layer <- sim_result$num_traps_per_layer

cat(sprintf("Number of layers: %d\n", n_layers))
cat(sprintf("Number of traps per layer: %s\n",
            paste(n_traps_per_layer, collapse=", ")))
cat(sprintf("Number of snapshots: %d\n\n", length(timepoints)))

# Summary statistics
cat("Total CO2 Volume Over Time:\n")
for (i in seq_along(timepoints)) {
  cat(sprintf("  Time %.1f years: %.2e m³ (%.2f Mt)\n",
              timepoints[i],
              total_volumes[i],
              total_volumes[i] * co2_density_l1 / 1e9))
}

cat("\n")

# Final distribution across layers
cat("Final CO2 Distribution by Layer:\n")
final_idx <- length(timepoints)
for (layer in 1:n_layers) {
  volume <- layer_volumes[final_idx, layer]
  percentage <- (volume / total_volumes[final_idx]) * 100
  cat(sprintf("  Layer %d: %.2e m³ (%.1f%%)\n", layer, volume, percentage))
}

cat("\n")

# ===========================
# Step 5: Visualize results
# ===========================
cat("=== Creating Visualizations ===\n\n")

# Create a simple plot of total CO2 volume over time
png("co2_total_volume.png", width=800, height=600)
plot(timepoints, total_volumes,
     type="b", pch=19, col="blue",
     xlab="Time (years)", ylab="Total CO2 Volume (m³)",
     main="Total CO2 Volume Over Time")
grid()
dev.off()
cat("  Saved: co2_total_volume.png\n")

# Create a stacked area plot showing layer contributions
png("co2_layer_distribution.png", width=800, height=600)
colors <- rainbow(n_layers)

# Prepare data for stacked plot
plot(timepoints, layer_volumes[, 1],
     type="n", ylim=c(0, max(total_volumes)),
     xlab="Time (years)", ylab="CO2 Volume (m³)",
     main="CO2 Distribution by Layer Over Time")

# Plot each layer as a polygon
bottom <- rep(0, length(timepoints))
for (layer in 1:n_layers) {
  top <- bottom + layer_volumes[, layer]
  polygon(c(timepoints, rev(timepoints)),
          c(top, rev(bottom)),
          col=colors[layer], border=NA)
  bottom <- top
}

# Add legend
legend("topleft", legend=paste("Layer", 1:n_layers),
       fill=colors, border=NA, cex=0.8)
grid()
dev.off()
cat("  Saved: co2_layer_distribution.png\n")

# Plot injection rate over time
png("injection_rate.png", width=800, height=600)
plot(injection_times, injection_rates_m3,
     type="h", lwd=3, col="darkgreen",
     xlab="Time (years)", ylab="Injection Rate (m³/year)",
     main="CO2 Injection Rate Over Time")
grid()
dev.off()
cat("  Saved: injection_rate.png\n\n")

cat("=== Example completed successfully! ===\n")
cat("\nThe simulation results show:\n")
cat(sprintf("  - Initial injection: %.2f Mt/year\n", annual_rates_mt[1]))
cat(sprintf("  - Peak injection: %.2f Mt/year\n", max(annual_rates_mt)))
cat(sprintf("  - Total injected: %.2f Mt over %d years\n",
            sum(annual_rates_mt), length(annual_rates_mt)))
cat(sprintf("  - Final storage: %.2f Mt in %d layers\n",
            total_volumes[final_idx] * co2_density_l1 / 1e9,
            n_layers))
cat(sprintf("  - Primary storage layer: Layer 1 (%.1f%% of total)\n",
            (layer_volumes[final_idx, 1] / total_volumes[final_idx]) * 100))

# ===========================
# Optional: Access trap-level data
# ===========================
cat("\n=== Trap-Level Analysis (Layer 1) ===\n\n")

# Get trap volumes for layer 1
trap_volumes_l1 <- sim_result$trap_co2_volumes[[1]]
trap_percentages_l1 <- sim_result$trap_co2_percentages[[1]]

# Find the traps with most CO2 at final time
n_traps_l1 <- n_traps_per_layer[1]
final_trap_volumes <- trap_volumes_l1[final_idx, ]
final_trap_percentages <- trap_percentages_l1[final_idx, ]

# Sort by volume
trap_order <- order(final_trap_volumes, decreasing=TRUE)
top_5_traps <- head(trap_order, 5)

cat("Top 5 traps by CO2 volume in Layer 1:\n")
for (i in seq_along(top_5_traps)) {
  trap_id <- top_5_traps[i]
  cat(sprintf("  Trap %d: %.2e m³ (%.1f%% of total)\n",
              trap_id,
              final_trap_volumes[trap_id],
              final_trap_percentages[trap_id]))
}

cat("\n=== All steps completed! ===\n")
