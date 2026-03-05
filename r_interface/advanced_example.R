 # CO2 Injection Simulation - Generic API Example
# Demonstrates the generic (dataset-agnostic) R interface for CO2BatchFill
# Uses setup_simulator_from_surfaces() to pass custom surface arrays directly

library(JuliaCall)

# Setup Julia and load package
julia_setup()

# REMEMBER TO SET WORKING DIR TO CO2BatchFill
# Alternatively: Install the package globally if you want to call the lib from
# another project.
julia_command('using Pkg; Pkg.activate(".")')

# This line might take a bit of time...
julia_command('using CO2BatchFill')

# ============================================================================
# ADVANCED EXAMPLE: Generic API with custom surface data
# ============================================================================

cat("=== Advanced Example: Generic API with Custom Surfaces ===\n\n")

# Step 1: Create synthetic surface data
# In practice, these would come from your own dataset (e.g., seismic interpretation)
cat("Creating synthetic surface data...\n")

nx <- 64
ny <- 64
dx <- 50.0   # 50m grid spacing
dy <- 50.0

# Create synthetic top and base surfaces for 3 layers
# Layer 1 (deepest): a gently dipping surface with some structural traps
x <- matrix(rep(1:nx, ny), nrow = nx, ncol = ny)
y <- matrix(rep(1:ny, each = nx), nrow = nx, ncol = ny)

# Base topography: gentle dip + dome structures
base_depth <- 1000.0 + 0.5 * x + 0.3 * y
dome1 <- 15.0 * exp(-((x - 20)^2 + (y - 30)^2) / 200)
dome2 <- 10.0 * exp(-((x - 45)^2 + (y - 40)^2) / 150)

layer1_top <- base_depth - dome1 - dome2  # top of layer 1 (structural traps at domes)
layer1_base <- base_depth + 20.0          # 20m thick sand

layer2_top <- base_depth - 30.0 - dome1 * 0.8 - dome2 * 0.5  # shallower layer
layer2_base <- layer1_top - 5.0            # 5m shale between layers

layer3_top <- base_depth - 60.0 - dome1 * 0.6   # shallowest layer
layer3_base <- layer2_top - 5.0                   # 5m shale between layers

# Step 2: Setup simulator using generic API
cat("Setting up simulator from surface arrays...\n")
setup_result <- julia_call("setup_simulator_from_surfaces",
                           layer_tops = list(layer1_top, layer2_top, layer3_top),
                           layer_bases = list(layer1_base, layer2_base, layer3_base),
                           layer_names = c("Deep", "Middle", "Shallow"),
                           dx = dx,
                           dy = dy,
                           boundary_condition = "closed")
print(setup_result)

nx_bc <- setup_result$nx_after_bc
ny_bc <- setup_result$ny_after_bc
n_layers <- setup_result$n_layers

# Step 3: Configure custom reservoir properties
cat("\nConfiguring custom reservoir properties...\n")

# Set shale pressure threshold — leakage height is computed automatically
# from shale_pressure_threshold / ((brine_density - co2_density) * g)
shale_threshold <- 80000.0  # Pa
# Make top layer impermeable (represents caprock)
shale_thresholds <- c(rep(shale_threshold, n_layers - 1), Inf)

config_result <- julia_call("configure_reservoir",
                             porosity = rep(0.35, n_layers),
                             residual_co2_sat = rep(0.15, n_layers),
                             irreducible_water_sat = rep(0.25, n_layers),
                             shale_pressure_threshold = shale_thresholds,
                             residual_leakage_time = rep(2.0, n_layers),
                             brine_density = rep(1020.0, n_layers),
                             co2_density = rep(500.0, n_layers),
                             layer_specific = TRUE)
print(config_result)

# Step 4: Setup injection scenario
cat("\nSetting up injection scenario...\n")

n_times <- 10
co2_density <- 500.0  # kg/m³

# Inject at dome1 location in layer 1
injection_rate_mt <- 0.5  # Mt/year constant
injection_rate_m3 <- injection_rate_mt * 1e9 / co2_density

layer1_injection <- array(0, dim = c(n_times, nx_bc, ny_bc))
for (i in 1:n_times) {
  # Inject near dome1 center (offset by 1 for closed BC padding)
  layer1_injection[i, 21, 31] <- injection_rate_m3
}

zero_injection <- array(0, dim = c(1, nx_bc, ny_bc))

injection_matrices <- list(
  layer1_injection,  # Layer 1 (deepest)
  zero_injection,    # Layer 2
  zero_injection     # Layer 3 (shallowest)
)

# Step 5: Run simulation
cat("\nRunning simulation...\n")
sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 10.0,
                         time_step = 1.0,
                         injection_rate_matrices = injection_matrices,
                         verbose = FALSE)

# Step 6: Check results
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

# Plot 1: Total CO2 volume over time
cat("\nPlotting total CO2 volume over time...\n")
timepoints <- sim_result$timepoints
total_volumes <- sim_result$total_co2_volumes
plot(timepoints, total_volumes / 1e6, type = "b",
     xlab = "Time (years)", ylab = "Total CO2 Volume (million m³)",
     main = "Total CO2 Volume Over Time (Generic API)",
     col = "blue", pch = 19)
grid()

# Plot 2: Layer-wise CO2 volumes
cat("\nPlotting layer-wise CO2 volumes...\n")
layer_volumes <- sim_result$layer_co2_volumes
matplot(timepoints, layer_volumes / 1e6, type = "b",
        xlab = "Time (years)", ylab = "CO2 Volume per Layer (million m³)",
        main = "Layer-wise CO2 Volumes (Generic API)",
        col = rainbow(n_layers), pch = 19, lty = 1)
legend("topleft", legend = c("Deep", "Middle", "Shallow"),
       col = rainbow(n_layers), pch = 19, lty = 1, cex = 0.8)
grid()

# Step 7: Generate animation
cat("\nGenerating animation...\n")
anim_result <- julia_call("generate_birdseye_animation",
                          output_file = "co2_generic_animation.gif",
                          num_frames = 20L,
                          fps = 2L,
                          max_co2_height = 15.0)

if (anim_result$status == "success") {
  cat("Animation saved to:", anim_result$output_file, "\n")
} else {
  cat("Animation failed:", anim_result$message, "\n")
  if (!is.null(anim_result$stacktrace)) {
    cat("Stacktrace:\n", anim_result$stacktrace, "\n")
  }
}
