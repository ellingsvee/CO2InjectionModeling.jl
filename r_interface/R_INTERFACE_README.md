# Using CO2InjectionModeling.jl from R

Use the CO2 injection simulator from R using the JuliaCall package.

## Setup

### 1. Install prerequisites

- Julia 1.6+ from [julialang.org](https://julialang.org/downloads/)
- R 4.0+ from [r-project.org](https://www.r-project.org/)
- JuliaCall R package:
  ```R
  install.packages("JuliaCall")
  ```

### 2. Install CO2InjectionModeling.jl

#### Option A: Use from the project directory

Navigate to the project directory and run:
```R
library(JuliaCall)
julia_setup()

# REMEMBER TO SET WORKING DIR TO CO2InjectionModeling.jl
julia_command('using Pkg; Pkg.activate(".")')

# This might take some time
julia_command('using CO2InjectionModeling')
```

#### Option B: Install as a package 

Install the package globally so you can use it from any directory:

```bash
cd /path/to/CO2InjectionModeling.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -e 'using Pkg; Pkg.develop(path=pwd())'
```

Then in R:
```R
library(JuliaCall)
julia_setup()
julia_command('using CO2InjectionModeling')
julia_command('using .CO2RInterface')
```


## Basic usage

```R
# 1. Setup simulator
julia_call("setup_simulator", boundary_condition = "open")

# 2. Configure reservoir
julia_call("configure_reservoir",
           porosity = 0.4,
           residual_co2_sat = 0.2,
           irreducible_water_sat = 0.3,
           shale_pressure_threshold = 98000.0,
           brine_co2_density_diff = 450.0,
           residual_leakage_time = 1.0,
           layer_specific = FALSE)

# 3. Run simulation
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

# 4. Access results
print(result$timepoints)
print(result$total_co2_volumes)
```

## Key parameters

### `setup_simulator()`
- `data_path`: Path to depth surfaces (default: `"sleipner/depth_surfaces/"`)
- `boundary_condition`: `"open"` or `"closed"` (default: `"open"`)

### `configure_reservoir()`
- `porosity`: Sand porosity (0-1)
- `residual_co2_sat`: Residual CO2 saturation (0-1)
- `irreducible_water_sat`: Irreducible water saturation (0-1)
- `shale_pressure_threshold`: Shale pressure threshold (Pa)
- `brine_co2_density_diff`: Density difference (kg/m³)
- `residual_leakage_time`: Residual leakage time (years)
- `layer_specific`: Set to `TRUE` to provide vectors for each layer

### `run_simulation()`
- `start_time`, `end_time`: Simulation time range (years)
- `time_step`: Time between snapshots (years)
- `injection_times`: Vector of injection times (years)
- `injection_locations_i`, `injection_locations_j`: Grid indices (1-based)
- `injection_amounts`: Injection rates (m³/year)
- `injection_layer_indices`: Layer indices (1-based, 1=bottom)
- `num_snapshots`: Number of snapshots to save
- `verbose`: Print progress

## Result structure

- `timepoints`: Vector of snapshot times
- `total_co2_volumes`: Total CO2 volume at each timepoint (m³)
- `layer_co2_volumes`: Matrix (timepoints × layers) of volumes per layer
- `trap_co2_volumes`: List of matrices for trap volumes
- `trap_co2_percentages`: List of matrices for trap percentages

## Examples

See example scripts:
- [r_interface_simple.R](r_interface_simple.R): Minimal example
- [r_interface_example.R](r_interface_example.R): Complete example with visualization
