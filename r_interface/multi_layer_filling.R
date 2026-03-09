# Multi-Layer CO2 Filling Example
# R equivalent of examples/multi_layer_filling.jl
#
# Creates synthetic 3-layer geology with dome-shaped structural traps,
# injects CO2 in the deepest layer, and tracks migration + leakage upward.

library(JuliaCall)

# Setup Julia and load package
julia_setup()

# Option A: activate from project directory
# setwd("/path/to/CO2BatchFill")
julia_command('using Pkg; Pkg.activate("..")')

# Option B: if CO2BatchFill is installed as a package, just load it directly
julia_command('using CO2BatchFill')

# Grid / domain settings
NX <- 100L
NY <- 100L
LENGTH_X <- 1000.0
LENGTH_Y <- 1000.0
DX <- LENGTH_X / NX
DY <- LENGTH_Y / NY
LAYER_THICK <- 10.0     # m (vertical thickness of each sand layer)
PAD_WIDTH <- 2L
N_LAYERS <- 3L
RESIDUAL_TRAPPING <- 0.4

# Generate synthetic surfaces
# In the Julia example, GaussianRandomFields are used. Here we create similar
# dome-shaped surfaces using deterministic functions for reproducibility.

cat("Generating synthetic surfaces...\n")

x <- matrix(rep(1:NX, NY), nrow = NX, ncol = NY) * DX
y <- matrix(rep(1:NY, each = NX), nrow = NX, ncol = NY) * DY

# Create dome-shaped structural traps at different positions
dome1 <- 5.0 * exp(-((x - 300)^2 + (y - 400)^2) / (2 * 150^2))
dome2 <- 4.0 * exp(-((x - 600)^2 + (y - 300)^2) / (2 * 120^2))
dome3 <- 3.5 * exp(-((x - 450)^2 + (y - 700)^2) / (2 * 100^2))
dome4 <- 3.0 * exp(-((x - 800)^2 + (y - 600)^2) / (2 * 130^2))

# Regional dip
dip <- 0.002 * x + 0.001 * y

# Surfaces at increasing depth (note: positive = deeper)
# Layer 3 (shallowest)
surf_L3 <- 850.0 + dip - dome1 * 0.8 - dome2 * 0.6 - dome3 * 0.7 - dome4 * 0.5
# Layer 2
surf_L2 <- 900.0 + dip - dome1 * 0.9 - dome2 * 0.7 - dome3 * 0.5 - dome4 * 0.8
# Layer 1 (deepest, injection layer)
surf_L1 <- 950.0 + dip - dome1 - dome2 - dome3 - dome4

# Setup simulator from surfaces
# Layers must be passed shallowest-first, deepest-last.
# The simulator internally reverses them so layers[1] = deepest (injection layer).

cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator_from_surfaces",
                           layer_tops = list(surf_L3, surf_L2, surf_L1),
                           layer_bases = list(surf_L3 + LAYER_THICK,
                                              surf_L2 + LAYER_THICK,
                                              surf_L1 + LAYER_THICK),
                           layer_names = c("Storage layer 3",
                                           "Storage layer 2",
                                           "Storage layer 1"),
                           dx = DX,
                           dy = DY,
                           boundary_condition = "closed",
                           pad_width = PAD_WIDTH)

cat("Setup result:\n")
print(setup_result)

nx_bc <- setup_result$nx_after_bc
ny_bc <- setup_result$ny_after_bc
n_layers <- setup_result$n_layers

# Configure reservoir properties
# Layer-specific: layers 1-2 have finite shale threshold, layer 3 (top) is sealed.

cat("\nConfiguring reservoir properties...\n")

residual_trapping <- 0.4
shale_threshold <- 15000.0  # Pa
shale_thresholds <- c(shale_threshold, shale_threshold, Inf)  # Layer 3 sealed (caprock)

config_result <- julia_call("configure_reservoir",
                             porosity = rep(0.3, n_layers),
                             residual_co2_sat = rep(residual_trapping, n_layers),
                             irreducible_water_sat = rep(0.1, n_layers),
                             shale_pressure_threshold = shale_thresholds,
                             residual_leakage_time = rep(5.0, n_layers),
                             layer_specific = TRUE)
cat("config result:\n")
print(config_result)

# Setup injection scenario
# Single central well in layer 1 (deepest) only, constant rate for 10 years.

cat("\nSetting up injection scenario...\n")

TOTAL_RATE <- 25000.0   # m^3/year
INJECTION_END <- 10.0
TIME_STEP <- 1.0

# Injection times: 0 to 9 years (rate active), then 10 years (rate = 0 = stop)
n_inject_times <- as.integer(INJECTION_END / TIME_STEP) + 1L  # 11 time steps

# Layer 1: inject at center
layer1_injection <- array(0, dim = c(n_inject_times, nx_bc, ny_bc))

# Center cell (accounting for boundary padding)
cx <- NX %/% 2L + PAD_WIDTH
cy <- NY %/% 2L + PAD_WIDTH

for (i in 1:(n_inject_times - 1L)) {
  layer1_injection[i, cx, cy] <- TOTAL_RATE
}
# Last time step: rate = 0 (stop injection)
# Already zero from array initialization

# Layers 2-3: no direct injection
zero_injection <- array(0, dim = c(1L, nx_bc, ny_bc))

injection_matrices <- list(
  layer1_injection,  # Layer 1 (deepest, injection layer)
  zero_injection,    # Layer 2 (CO2 arrives via leakage from below)
  zero_injection     # Layer 3 (CO2 arrives via leakage from below)
)

# Run simulation
cat("\nRunning multi-layer fill simulation...\n")

sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 15.0,
                         time_step = TIME_STEP,
                         injection_rate_matrices = injection_matrices,
                         verbose = FALSE)

if (sim_result$status == "success") {
  cat("\nSimulation Successful!\n")
  cat("Timepoints:", sim_result$timepoints, "\n")
  cat("Number of layers:", sim_result$num_layers, "\n")
  cat("Traps per layer:", sim_result$num_traps_per_layer, "\n")
  cat("\nFinal total CO2 volume:", tail(sim_result$total_co2_volumes, 1), "m^3\n")
} else {
  cat("\nSimulation failed:", sim_result$message, "\n")
  if (!is.null(sim_result$stacktrace)) {
    cat("Stacktrace:\n", sim_result$stacktrace, "\n")
  }
  stop("Simulation failed")
}

# Plot results
timepoints <- sim_result$timepoints
total_volumes <- sim_result$total_co2_volumes
layer_volumes <- sim_result$layer_co2_volumes

# Plot 1: Total CO2 volume over time
cat("\nPlotting total CO2 volume over time...\n")
plot(timepoints, total_volumes,
     type = "b", pch = 19, col = "blue",
     xlab = "Time (years)",
     ylab = expression("Total CO2 Volume"),
     main = "Total CO2 Volume Over Time")
grid()

# Plot 2: CO2 volume per layer over time
cat("Plotting layer-wise CO2 volumes...\n")
colors <- c("blue", "red", "darkgreen")
layer_labels <- c("Storage layer 1", "Storage layer 2", "Storage layer 3")

matplot(timepoints, layer_volumes,
        type = "b", pch = 19, lty = 1,
        col = colors,
        xlab = "Time (years)",
        ylab = expression("CO2 Volume per Layer"),
        main = "Layer-wise CO2 Volumes Over Time")
legend("topleft", legend = layer_labels,
       col = colors, pch = 19, lty = 1, cex = 0.8)
grid()