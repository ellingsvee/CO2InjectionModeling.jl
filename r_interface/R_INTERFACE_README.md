# Using CO2BatchFill.jl from R

Use the CO2 injection simulator from R using the JuliaCall package.

## Setup

### 1. Install prerequisites

- Julia 1.6+ from [julialang.org](https://julialang.org/downloads/)
- R 4.0+ from [r-project.org](https://www.r-project.org/)
- JuliaCall R package:
  ```R
  install.packages("JuliaCall")
  ```

### 2. Install CO2BatchFill.jl

#### Option A: Use from the project directory

Navigate to the project directory and run:
```R
library(JuliaCall)
julia_setup()

# REMEMBER TO SET WORKING DIR TO CO2BatchFill
julia_command('using Pkg; Pkg.activate(".")')

# This might take some time
julia_command('using CO2BatchFill')
```

#### Option B: Install as a package

Install the package globally so you can use it from any directory:

```bash
cd /path/to/CO2BatchFill
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -e 'using Pkg; Pkg.develop(path=pwd())'
```

Then in R:
```R
library(JuliaCall)
julia_setup()
julia_command('using CO2BatchFill')
```

## Two APIs

CO2BatchFill provides two ways to run simulations from R:

### Generic API (any dataset)

Use `setup_simulator_from_surfaces()` to pass your own surface arrays directly. No dataset-specific code needed.

```R
julia_command('using CO2BatchFill')

setup <- julia_call("setup_simulator_from_surfaces",
                    layer_tops = list(top_surface_1, top_surface_2),
                    layer_bases = list(base_surface_1, base_surface_2),
                    layer_names = c("Layer1", "Layer2"),
                    dx = 50.0, dy = 50.0,
                    boundary_condition = "closed")
```

### Sleipner Convenience API

For Sleipner data, use the convenience wrappers that load the data and set default reservoir properties:

```R
julia_command('using CO2BatchFill')
julia_command('include("examples/Sleipner/Sleipner.jl"); using .Sleipner; using .Sleipner.CO2RInterface')

setup <- julia_call("setup_sleipner_simulator", boundary_condition = "closed")
julia_call("setup_sleipner_reservoir")
```

## Quick Start (Generic API)

```R
library(JuliaCall)
julia_setup()
julia_command('using Pkg; Pkg.activate(".")')
julia_command('using CO2BatchFill')

# 1. Create or load your surface data (nx × ny matrices)
nx <- 64; ny <- 64
top1 <- matrix(runif(nx*ny, 900, 950), nrow=nx, ncol=ny)
base1 <- top1 + 20  # 20m thick layer

# 2. Setup simulator
setup <- julia_call("setup_simulator_from_surfaces",
                    layer_tops = list(top1),
                    layer_bases = list(base1),
                    layer_names = c("L1"),
                    dx = 50.0, dy = 50.0,
                    boundary_condition = "closed")

nx_bc <- setup$nx_after_bc
ny_bc <- setup$ny_after_bc

# 3. Configure reservoir properties
julia_call("configure_reservoir",
           porosity = 0.35,
           residual_co2_sat = 0.2,
           irreducible_water_sat = 0.3,
           shale_pressure_threshold = 98000.0,
           leakage_height = Inf,
           residual_leakage_time = 1.0)

# 4. Create injection scenario
n_times <- 10
injection <- array(0, dim = c(n_times, nx_bc, ny_bc))
for (i in 1:n_times) {
  injection[i, 32, 32] <- 1e6  # m³/year at center
}
injection_matrices <- list(injection)

# 5. Run simulation
result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 10.0,
                     time_step = 1.0,
                     injection_rate_matrices = injection_matrices)

print(result$total_co2_volumes)
```

## Quick Start (Sleipner)

```R
library(JuliaCall)
julia_setup()
julia_command('using Pkg; Pkg.activate(".")')
julia_command('using CO2BatchFill')
julia_command('include("examples/Sleipner/Sleipner.jl"); using .Sleipner; using .Sleipner.CO2RInterface')

# 1. Setup simulator with Sleipner data
setup <- julia_call("setup_sleipner_simulator", boundary_condition = "open")
nx <- setup$nx
ny <- setup$ny

# 2. Use Sleipner default properties
julia_call("setup_sleipner_reservoir")

# 3. Create injection scenario
n_times <- 15
layer1_injection <- array(0, dim = c(n_times, nx, ny))

rates_mt <- c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
              0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
              0.82, 0.86, 0.76)  # Mt/year
rates_m3 <- rates_mt * 1e9 / 570  # Convert to m³/year

for (i in 1:n_times) {
  layer1_injection[i, 32, 59] <- rates_m3[i]
}

zero_injection <- array(0, dim = c(1, nx, ny))

injection_matrices <- list(
  layer1_injection,
  zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection
)

# 4. Run simulation
result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 15.0,
                     time_step = 1.0,
                     injection_rate_matrices = injection_matrices,
                     verbose = TRUE)

print(result$timepoints)
print(result$total_co2_volumes)
```

## API Reference

### Generic API (CO2BatchFill)

#### `setup_simulator_from_surfaces()`

Setup the simulator from raw surface arrays. This is the primary entry point for custom datasets.

**Parameters:**
- `layer_tops`: List of 2D matrices (each nx × ny), top surfaces for each layer
- `layer_bases`: List of 2D matrices (each nx × ny), base surfaces for each layer
- `layer_names`: Character vector of layer names (e.g., `c("L1", "L2")`)
- `dx`: Grid spacing in x direction (meters)
- `dy`: Grid spacing in y direction (meters)
- `boundary_condition`: `"open"` or `"closed"` (default: `"open"`)
- `caprock_surface`: Optional caprock top surface matrix (nx × ny), or `NULL`

**Returns:** Dictionary with `status`, `n_layers`, `nx`, `ny`, `boundary_condition`, `nx_after_bc`, `ny_after_bc`

---

#### `configure_reservoir()`

Configure custom reservoir properties.

**Parameters:**
- `porosity`: Sand porosity (0-1). Scalar or vector of length n_layers
- `residual_co2_sat`: Residual CO2 saturation (0-1). Scalar or vector
- `irreducible_water_sat`: Irreducible water saturation (0-1). Scalar or vector
- `shale_pressure_threshold`: Shale pressure threshold (Pa). Scalar or vector
- `leakage_height`: Critical CO2 height for leakage through shale (m). Scalar or vector. Use `Inf` for impermeable caprock
- `residual_leakage_time`: Residual leakage time (years). Scalar or vector
- `layer_specific`: Set to `TRUE` to provide vectors for layer-specific properties (default: `FALSE`)

**Returns:** Dictionary with `status`, `n_layers`

**Example (uniform properties):**
```R
julia_call("configure_reservoir",
           porosity = 0.35,
           residual_co2_sat = 0.2,
           irreducible_water_sat = 0.3,
           shale_pressure_threshold = 98000.0,
           leakage_height = 17.0,
           residual_leakage_time = 1.0,
           layer_specific = FALSE)
```

**Example (layer-specific properties):**
```R
julia_call("configure_reservoir",
           porosity = rep(0.4, 3),
           residual_co2_sat = rep(0.2, 3),
           irreducible_water_sat = rep(0.3, 3),
           shale_pressure_threshold = rep(98000.0, 3),
           leakage_height = c(17.0, 17.0, Inf),
           residual_leakage_time = rep(1.0, 3),
           layer_specific = TRUE)
```

---

#### `run_simulation()`

Run a CO2 injection simulation.

**Parameters:**
- `start_time`: Simulation start time (years)
- `end_time`: Simulation end time (years)
- `time_step`: Time step for output snapshots (years)
- `injection_rate_matrices`: List of 3D arrays (one per layer), each with dimensions `(n_times × nx × ny)`. Each array specifies injection rates (m³/year) at each grid cell for each time point. For layers without injection, provide a `(1 × nx × ny)` array of zeros.
- `verbose`: Print progress messages (default: `FALSE`)

**Returns:** Dictionary containing:
- `status`: "success" or "error"
- `timepoints`: Vector of snapshot times (years)
- `total_co2_volumes`: Total CO2 stored in the reservoir at each timepoint (m³)
- `layer_co2_volumes`: Matrix (timepoints × layers) of volumes per layer (m³)
- `num_layers`: Number of layers
- `num_traps_per_layer`: Vector of trap counts per layer

---

#### `generate_birdseye_animation()`

Generate a bird's eye view animation showing CO2 distribution across all layers over time.

**Parameters:**
- `output_file`: Path where to save the animation (default: `"multi_layer_filling.gif"`)
- `num_frames`: Number of frames in animation (default: `30L`)
- `start_time`: Start time for animation in years (default: `0.0`)
- `end_time`: End time for animation in years, or `NULL` for auto-detect (default: `NULL`)
- `fps`: Frames per second (default: `2L`)
- `colormap`: Colormap name for CO2 heights (default: `"thermal"`)
- `max_CO2_height`: Maximum CO2 height for colorscale in meters (default: `20.0`)

**Returns:** Dictionary with `status`, `output_file`

---

### Sleipner Convenience API

These functions are available after loading the Sleipner module:
```R
julia_command('include("examples/Sleipner/Sleipner.jl"); using .Sleipner; using .Sleipner.CO2RInterface')
```

#### `setup_sleipner_simulator()`

Load Sleipner topography and set up the simulator with Sleipner defaults.

**Parameters:**
- `data_path`: Path to depth surfaces (default: `"sleipner/depth_surfaces/"`)
- `boundary_condition`: `"open"` or `"closed"` (default: `"open"`)

**Returns:** Dictionary with `status`, `n_layers`, `nx`, `ny`, `boundary_condition`, `nx_after_bc`, `ny_after_bc`

---

#### `setup_sleipner_reservoir()`

Configure reservoir properties using Sleipner field default values:
- Porosity: 0.4
- Residual CO2 saturation: 0.2
- Irreducible water saturation: 0.3
- Shale pressure threshold: 98000.0 Pa
- Leakage heights: Computed from density differences (brine: 1020 kg/m³, CO2: 460 kg/m³)
- Residual leakage time: 5.0 years

**Parameters:** None (must call `setup_sleipner_simulator()` first)

**Returns:** Dictionary with `status`, `n_layers`

---

## Result Structure

All functions return a dictionary with a `status` field:
- `"success"`: Operation completed successfully
- `"error"`: Operation failed (check `message` field for details)

Simulation results include:
- `timepoints`: Times when snapshots were taken
- `total_co2_volumes`: Total CO2 stored in the reservoir at each timepoint
- `layer_co2_volumes`: CO2 volume breakdown by layer (timepoints × layers matrix)
- `num_layers`: Number of layers in the simulation
- `num_traps_per_layer`: Number of traps in each layer

## Examples

- [example.R](example.R) — Sleipner convenience API (loads Sleipner data, uses defaults)
- [advanced_example.R](advanced_example.R) — Generic API with custom surface arrays

## Notes

- **Grid indexing**: R uses 1-based indexing, which matches the Julia interface
- **Array dimensions**: Injection matrices use dimension order `(n_times × nx × ny)`
- **Units**:
  - Injection rates: m³/year
  - Volumes: m³
  - Time: years
  - Pressure: Pa
  - Density: kg/m³
- **Performance**: First simulation run will be slower due to Julia JIT compilation. Subsequent runs are much faster.
