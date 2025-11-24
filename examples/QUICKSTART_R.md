# Quick Start Guide: Using CO2InjectionModeling from R

## Installation

### 1. Install Prerequisites

```R
# Install JuliaCall
install.packages("JuliaCall")
```

### 2. Setup Julia Environment

From the CO2InjectionModeling.jl directory:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Your First Simulation

Create a new R script with the following code:

```R
library(JuliaCall)

# Initialize Julia (adjust path if needed)
julia_setup()

# Load the package
julia_command('using Pkg; Pkg.activate(".")')
julia_command('using CO2InjectionModeling')

# Step 1: Load Sleipner topography
cat("Loading Sleipner data...\n")
setup <- julia_call("setup_simulator",
                    boundary_condition = "open")
print(setup)

# Step 2: Configure reservoir properties
cat("\nConfiguring reservoir...\n")
config <- julia_call("configure_reservoir",
                     porosity = 0.4,
                     residual_co2_sat = 0.2,
                     irreducible_water_sat = 0.3,
                     shale_pressure_threshold = 98000.0,
                     brine_co2_density_diff = 450.0,
                     residual_leakage_time = 1.0,
                     layer_specific = FALSE)
print(config)

# Step 3: Run a simple simulation
cat("\nRunning simulation...\n")
result <- julia_call("run_simulation",
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

# Display results
cat("\n=== Results ===\n")
cat("Timepoints:", result$timepoints, "\n")
cat("Total CO2 volumes:", result$total_co2_volumes, "\n")

# Create a simple plot
plot(result$timepoints, result$total_co2_volumes,
     type = "b", pch = 19, col = "blue",
     xlab = "Time (years)",
     ylab = "CO2 Volume (m³)",
     main = "CO2 Storage Over Time")
grid()
```

## Running the Examples

The package includes ready-to-run examples:

### Simple Example
```bash
cd examples
Rscript r_interface_simple.R
```

### Comprehensive Example (with visualization)
```bash
cd examples
Rscript r_interface_example.R
```

This will create three plots:
- `co2_total_volume.png` - Total CO2 over time
- `co2_layer_distribution.png` - Distribution across layers
- `injection_rate.png` - Injection rate schedule

## Common Issues

### Julia not found
```R
julia_setup(JULIA_HOME = "/path/to/julia/bin")
```

### Package not loading
Make sure you're in the correct directory and the package is installed:
```R
julia_command('using Pkg; Pkg.activate("/full/path/to/CO2InjectionModeling.jl")')
```

### Memory issues
For large simulations:
```R
julia_setup(JULIA_FLAGS = "--heap-size-hint=4G")
```

## Next Steps

- Read [R_INTERFACE_README.md](R_INTERFACE_README.md) for full API documentation
- Modify the examples for your use case
- See [r_interface_example.R](r_interface_example.R) for advanced features

## Getting Help

- Check the [R_INTERFACE_README.md](R_INTERFACE_README.md) for detailed documentation
- Run the test suite: `julia --project=. examples/test_r_interface.jl`
- Look at the comprehensive example: [r_interface_example.R](r_interface_example.R)

## Understanding the Output

The simulation returns a dictionary with:

```R
result$timepoints              # Vector of times (years)
result$total_co2_volumes       # Total CO2 at each time (m³)
result$layer_co2_volumes       # Matrix: timepoints × layers
result$trap_co2_volumes        # List of matrices per layer
result$trap_co2_percentages    # Percentage of total in each trap
result$num_layers              # Number of layers (9 for Sleipner)
result$num_traps_per_layer     # Vector of trap counts
```

### Accessing Data

```R
# Get volume in layer 1 over time
layer1 <- result$layer_co2_volumes[, 1]

# Get final state (all layers)
final_idx <- length(result$timepoints)
final_volumes <- result$layer_co2_volumes[final_idx, ]

# Get trap data for layer 1
trap_volumes_l1 <- result$trap_co2_volumes[[1]]
```

## Example: Historical Sleipner Injection

```R
# Replicate 1996-2010 Sleipner injection
annual_rates_mt <- c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
                     0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
                     0.82, 0.86, 0.76)

# Convert Mt/year to m³/year
co2_density <- 570.0  # kg/m³ at bottom layer
injection_rates <- annual_rates_mt * 1e9 / co2_density

# Run simulation
result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 15.0,
                     time_step = 1.0,
                     injection_times = seq(0, 14),
                     injection_locations_i = rep(32L, 15),
                     injection_locations_j = rep(59L, 15),
                     injection_amounts = injection_rates,
                     injection_layer_indices = rep(1L, 15),
                     num_snapshots = 14L,
                     verbose = TRUE)
```

## Tips

1. **Start small** - Test with short simulations first
2. **Check status** - All functions return a `status` field
3. **Use verbose** - Set `verbose = TRUE` to see progress
4. **Save results** - Store results as RDS files for later analysis
   ```R
   saveRDS(result, "simulation_results.rds")
   ```
5. **Visualize** - R has excellent plotting capabilities - use them!

## Performance Notes

- First simulation is slower (Julia compilation)
- Subsequent runs are much faster
- State is maintained between calls
- Large grids or long time periods will take longer
- Consider reducing `num_snapshots` for faster runs

## Converting Units

Common conversions:

```R
# Mt to m³ (using density)
volume_m3 <- mass_mt * 1e9 / density_kg_m3

# m³ to Mt
mass_mt <- volume_m3 * density_kg_m3 / 1e9

# Typical CO2 densities (kg/m³) at Sleipner depths:
# L1 (bottom): 570
# L9 (top): 350
```

Happy simulating! 🎉
