# R Interface for CO2InjectionModeling.jl

This document describes how to use the CO2 injection simulator from R using the JuliaCall package.

## Installation

### Prerequisites

1. **Julia**: Install Julia (version 1.6+) from [julialang.org](https://julialang.org/downloads/)
2. **R**: Install R (version 4.0+) from [r-project.org](https://www.r-project.org/)
3. **JuliaCall**: Install the R package:
   ```R
   install.packages("JuliaCall")
   ```

### Setup

1. Clone and setup the CO2InjectionModeling.jl package
2. Make sure the Sleipner data is available in `sleipner/depth_surfaces/`

## Quick Start

```R
library(JuliaCall)

# Initialize Julia
julia_setup()
julia_command('using Pkg; Pkg.activate(".")')
julia_command('using CO2InjectionModeling')
julia_command('using .CO2RInterface')

# Step 1: Setup simulator
setup_result <- julia_call("setup_simulator",
                           boundary_condition = "open")

# Step 2: Configure reservoir properties
config_result <- julia_call("configure_reservoir",
                            porosity = 0.4,
                            residual_co2_sat = 0.2,
                            irreducible_water_sat = 0.3,
                            shale_pressure_threshold = 98000.0,
                            brine_co2_density_diff = 450.0,
                            residual_leakage_time = 1.0,
                            layer_specific = FALSE)

# Step 3: Run simulation
sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 5.0,
                         time_step = 1.0,
                         injection_times = c(0, 1, 2, 3, 4),
                         injection_locations_i = rep(32L, 5),
                         injection_locations_j = rep(59L, 5),
                         injection_amounts = rep(1e6, 5),
                         injection_layer_indices = rep(1L, 5),
                         num_snapshots = 5L,
                         verbose = TRUE)

# Access results
print(sim_result$timepoints)
print(sim_result$total_co2_volumes)
```

## API Reference

### 1. `setup_simulator()`

Loads and sets up the Sleipner topography layers.

**Parameters:**
- `data_path` (string, optional): Path to depth surfaces data. Default: `"sleipner/depth_surfaces/"`
- `boundary_condition` (string, optional): Either `"open"` or `"closed"`. Default: `"open"`

**Returns:** Dictionary with:
- `status`: `"success"` or `"error"`
- `n_layers`: Number of layers loaded
- `nx`, `ny`: Grid dimensions
- `boundary_condition`: The boundary condition used

**Example:**
```R
result <- julia_call("setup_simulator",
                     data_path = "sleipner/depth_surfaces/",
                     boundary_condition = "open")
```

### 2. `configure_reservoir()`

Configures reservoir properties for all layers.

**Parameters:**
- `porosity` (float or vector): Sand porosity (0-1)
- `residual_co2_sat` (float or vector): Residual CO2 saturation (0-1)
- `irreducible_water_sat` (float or vector): Irreducible water saturation (0-1)
- `shale_pressure_threshold` (float or vector): Shale pressure threshold (Pa)
- `brine_co2_density_diff` (float or vector): Density difference between brine and CO2 (kg/m³)
- `residual_leakage_time` (float or vector): Residual leakage time (years)
- `layer_specific` (boolean, optional): If TRUE, expects vectors of length n_layers. Default: FALSE

**Returns:** Dictionary with:
- `status`: `"success"` or `"error"`
- `n_layers`: Number of layers configured

**Example (uniform properties):**
```R
result <- julia_call("configure_reservoir",
                     porosity = 0.4,
                     residual_co2_sat = 0.2,
                     irreducible_water_sat = 0.3,
                     shale_pressure_threshold = 98000.0,
                     brine_co2_density_diff = 450.0,
                     residual_leakage_time = 1.0,
                     layer_specific = FALSE)
```

**Example (layer-specific properties):**
```R
# Different properties for each of the 9 Sleipner layers
result <- julia_call("configure_reservoir",
                     porosity = rep(0.4, 9),
                     residual_co2_sat = rep(0.2, 9),
                     irreducible_water_sat = rep(0.3, 9),
                     shale_pressure_threshold = rep(98000.0, 9),
                     brine_co2_density_diff = c(450, 477.5, 505, 532.5, 560,
                                                587.5, 615, 642.5, 670),
                     residual_leakage_time = rep(1.0, 9),
                     layer_specific = TRUE)
```

### 3. `run_simulation()`

Runs the CO2 injection simulation.

**Parameters:**
- `start_time` (float): Simulation start time (years)
- `end_time` (float): Simulation end time (years)
- `time_step` (float): Time step between snapshots (years)
- `injection_times` (vector of floats): Times of injection events (years)
- `injection_locations_i` (vector of integers): i-indices (1-based) for injection locations
- `injection_locations_j` (vector of integers): j-indices (1-based) for injection locations
- `injection_amounts` (vector of floats): Injection rates (m³/year)
- `injection_layer_indices` (vector of integers): Layer indices (1-based, 1=bottom layer)
- `num_snapshots` (integer, optional): Number of snapshots to save
- `verbose` (boolean, optional): Print progress messages. Default: FALSE

**Returns:** Dictionary with:
- `status`: `"success"` or `"error"`
- `timepoints`: Vector of snapshot times
- `total_co2_volumes`: Total CO2 volume at each timepoint (m³)
- `layer_co2_volumes`: Matrix (timepoints × layers) of CO2 volumes per layer (m³)
- `trap_co2_volumes`: List of matrices for trap volumes in each layer
- `trap_co2_percentages`: List of matrices for trap percentages
- `num_layers`: Number of layers
- `num_traps_per_layer`: Vector of trap counts per layer

**Example:**
```R
result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 15.0,
                     time_step = 1.0,
                     injection_times = seq(0, 14),
                     injection_locations_i = rep(32L, 15),
                     injection_locations_j = rep(59L, 15),
                     injection_amounts = rep(1e6, 15),
                     injection_layer_indices = rep(1L, 15),
                     num_snapshots = 14L,
                     verbose = TRUE)
```

## Understanding the Coordinate System

- **Layers**: Numbered from 1 (bottom) to 9 (top) for Sleipner
- **Grid indices**: 1-based indexing (Julia convention accessible from R)
  - `i`: Row index (1 to nx)
  - `j`: Column index (1 to ny)
- **Time**: In years
- **Volumes**: In cubic meters (m³)
- **Injection rates**: In m³/year

## Working with Results

### Accessing Data

```R
# After running simulation
sim_result <- julia_call("run_simulation", ...)

# Time series data
times <- sim_result$timepoints
total_vols <- sim_result$total_co2_volumes

# Layer-specific data (matrix: timepoints × layers)
layer_vols <- sim_result$layer_co2_volumes

# Get volume in layer 1 at all timepoints
layer1_vols <- layer_vols[, 1]

# Get volumes in all layers at final timepoint
final_idx <- length(times)
final_layer_vols <- layer_vols[final_idx, ]

# Trap-level data (list of matrices)
trap_vols_layer1 <- sim_result$trap_co2_volumes[[1]]  # Timepoints × traps
trap_pcts_layer1 <- sim_result$trap_co2_percentages[[1]]
```

### Plotting Results

```R
# Plot total volume over time
plot(sim_result$timepoints, sim_result$total_co2_volumes,
     type="b", xlab="Time (years)", ylab="CO2 Volume (m³)")

# Plot layer contributions
matplot(sim_result$timepoints, sim_result$layer_co2_volumes,
        type="l", xlab="Time (years)", ylab="CO2 Volume (m³)",
        main="CO2 Distribution by Layer")
legend("topleft", legend=paste("Layer", 1:sim_result$num_layers),
       col=1:sim_result$num_layers, lty=1:sim_result$num_layers)
```

## Example Use Cases

### Case 1: Historical Sleipner Injection

Replicate the historical Sleipner injection (1996-2010):

```R
# Historical injection rates (Mt/year)
annual_rates_mt <- c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
                     0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
                     0.82, 0.86, 0.76)

# Convert to m³/year (CO2 density = 570 kg/m³)
injection_rates <- annual_rates_mt * 1e9 / 570

result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 15.0,
                     time_step = 1.0,
                     injection_times = seq(0, 14),
                     injection_locations_i = rep(32L, 15),
                     injection_locations_j = rep(59L, 15),
                     injection_amounts = injection_rates,
                     injection_layer_indices = rep(1L, 15),
                     num_snapshots = 14L)
```

### Case 2: Multiple Injection Wells

Inject at different locations:

```R
# Two wells injecting simultaneously
injection_times <- c(0, 0, 1, 1, 2, 2)  # Each time repeated for each well
injection_i <- c(30L, 40L, 30L, 40L, 30L, 40L)  # Alternate between wells
injection_j <- c(50L, 50L, 50L, 50L, 50L, 50L)
injection_amounts <- c(5e5, 5e5, 5e5, 5e5, 5e5, 5e5)  # 0.5 million m³/year each
injection_layers <- rep(1L, 6)  # All in bottom layer

result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 3.0,
                     time_step = 1.0,
                     injection_times = injection_times,
                     injection_locations_i = injection_i,
                     injection_locations_j = injection_j,
                     injection_amounts = injection_amounts,
                     injection_layer_indices = injection_layers,
                     num_snapshots = 3L)
```

### Case 3: Sensitivity Analysis

Run multiple simulations with different parameters:

```R
# Test different porosities
porosities <- seq(0.3, 0.5, by=0.05)
results <- list()

for (i in seq_along(porosities)) {
  # Setup and configure
  julia_call("setup_simulator", boundary_condition = "open")
  julia_call("configure_reservoir",
             porosity = porosities[i],
             residual_co2_sat = 0.2,
             irreducible_water_sat = 0.3,
             shale_pressure_threshold = 98000.0,
             brine_co2_density_diff = 450.0,
             residual_leakage_time = 1.0)

  # Run simulation
  results[[i]] <- julia_call("run_simulation",
                             start_time = 0.0,
                             end_time = 10.0,
                             time_step = 1.0,
                             injection_times = seq(0, 9),
                             injection_locations_i = rep(32L, 10),
                             injection_locations_j = rep(59L, 10),
                             injection_amounts = rep(1e6, 10),
                             injection_layer_indices = rep(1L, 10),
                             num_snapshots = 10L)
}

# Compare results
for (i in seq_along(porosities)) {
  cat(sprintf("Porosity %.2f: Final volume = %.2e m³\n",
              porosities[i],
              tail(results[[i]]$total_co2_volumes, 1)))
}
```

## Troubleshooting

### Julia not found
```R
# Specify Julia path explicitly
julia_setup(JULIA_HOME = "/path/to/julia/bin")
```

### Package not found
```R
# Make sure to activate the project
julia_command('using Pkg; Pkg.activate("/path/to/CO2InjectionModeling.jl")')
julia_command('Pkg.instantiate()')
```

### Memory issues
For large simulations, you may need to increase Julia's memory:
```R
julia_setup(JULIA_HOME = "/path/to/julia/bin",
            JULIA_FLAGS = "--heap-size-hint=4G")
```

## Notes

- All indices are 1-based (Julia/R convention)
- Time is always in years
- Volumes are in cubic meters (m³)
- Injection rates are in m³/year
- The simulator state is maintained between calls, so you only need to call `setup_simulator()` and `configure_reservoir()` once for multiple simulations with the same configuration
- To change configuration, call `setup_simulator()` again to reset

## Examples

See the example scripts:
- `r_interface_simple.R`: Minimal working example
- `r_interface_example.R`: Complete example with visualization and analysis
