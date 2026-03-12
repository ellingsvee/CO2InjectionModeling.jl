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

# Set working directory to CO2BatchFill root
julia_command('using Pkg; Pkg.activate(".")')

# This might take some time on first run
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

## Quick Start

```R
library(JuliaCall)
julia_setup()
julia_command('using CO2BatchFill')

# 1. Create or load surface data (nx x ny matrices, shallowest-first)
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
           shale_pressure_threshold = Inf,
           residual_leakage_time = 1.0)

# 4. Create injection scenario
n_times <- 10
injection <- array(0, dim = c(n_times, nx_bc, ny_bc))
for (i in 1:n_times) {
  injection[i, 32 + 1, 32 + 1] <- 1e6  # m^3/year at center (offset by pad_width = 1)
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

## API Reference

### `setup_simulator_from_surfaces()`

Setup the simulator from raw surface arrays. This is the primary entry point for custom datasets.

**Parameters:**
- `layer_tops`: List of 2D matrices (each nx x ny), top surfaces for each layer. **Ordered shallowest-first, deepest-last.**
- `layer_bases`: List of 2D matrices (each nx x ny), base surfaces for each layer
- `layer_names`: Character vector of layer names (e.g., `c("L1", "L2")`)
- `dx`: Grid spacing in x direction (meters)
- `dy`: Grid spacing in y direction (meters)
- `boundary_condition`: `"open"` or `"closed"` (default: `"open"`)
- `caprock_surface`: Optional caprock top surface matrix (nx x ny), or `NULL`

**Returns:** Dictionary with `status`, `n_layers`, `nx`, `ny`, `boundary_condition`, `nx_after_bc`, `ny_after_bc`

**Note on layer ordering:** Layers must be provided shallowest-first, deepest-last. The simulator internally reverses them so `layers[1]` = deepest (injection layer). Leakage propagates from deeper layers upward.

---

### `configure_reservoir()`

Configure reservoir properties.

**Parameters:**
- `porosity`: Sand porosity (0-1). Scalar or vector of length n_layers
- `residual_co2_sat`: Residual CO2 saturation (0-1). Scalar or vector
- `irreducible_water_sat`: Irreducible water saturation (0-1). Scalar or vector
- `shale_pressure_threshold`: Shale pressure threshold (Pa). Scalar or vector. Use `Inf` for impermeable caprock
- `residual_leakage_time`: Residual leakage time (years). Scalar or vector
- `brine_density`: Brine density (kg/m^3). Scalar or vector (default: 1020.0)
- `co2_density`: CO2 density (kg/m^3). Scalar or vector (default: 460.0)
- `layer_specific`: Set to `TRUE` to provide vectors for layer-specific properties (default: `FALSE`)

**Returns:** Dictionary with `status`, `n_layers`

**Example (uniform properties):**
```R
julia_call("configure_reservoir",
           porosity = 0.35,
           residual_co2_sat = 0.2,
           irreducible_water_sat = 0.3,
           shale_pressure_threshold = 98000.0,
           residual_leakage_time = 1.0)
```

**Example (layer-specific, sealed top layer):**
```R
julia_call("configure_reservoir",
           porosity = rep(0.3, 3),
           residual_co2_sat = rep(0.4, 3),
           irreducible_water_sat = rep(0.1, 3),
           shale_pressure_threshold = c(15000.0, 15000.0, Inf),
           residual_leakage_time = rep(5.0, 3),
           layer_specific = TRUE)
```

---

### `run_simulation()`

Run a CO2 injection simulation.

**Parameters:**
- `start_time`: Simulation start time (years)
- `end_time`: Simulation end time (years)
- `time_step`: Time step for output snapshots (years)
- `injection_rate_matrices`: List of 3D arrays (one per layer), each with dimensions `(n_times x nx_after_bc x ny_after_bc)`. Each array specifies injection rates (m^3/year) at each grid cell for each time point. For layers without injection, provide a `(1 x nx_after_bc x ny_after_bc)` array of zeros.
- `verbose`: Print progress messages (default: `FALSE`)

**Returns:** Dictionary containing:
- `status`: "success" or "error"
- `timepoints`: Vector of snapshot times (years)
- `total_co2_volumes`: Total CO2 stored in the reservoir at each timepoint (m^3)
- `layer_co2_volumes`: Matrix (timepoints x layers) of volumes per layer (m^3)
- `num_layers`: Number of layers
- `num_traps_per_layer`: Vector of trap counts per layer

---

### `generate_birdseye_animation()`

Generate a bird's eye view animation showing CO2 distribution across all layers over time. Requires CairoMakie to be available in the Julia environment.

**Parameters:**
- `output_file`: Path where to save the animation (default: `"multi_layer_filling.gif"`)
- `num_frames`: Number of frames in animation (default: `30L`)
- `start_time`: Start time for animation in years (default: `NULL` for auto)
- `end_time`: End time for animation in years (default: `NULL` for auto)
- `fps`: Frames per second (default: `2L`)
- `colormap`: Colormap name for CO2 heights (default: `"thermal"`)
- `max_co2_height`: Maximum CO2 height for colorscale in meters (default: `20.0`)

**Returns:** Dictionary with `status`, `output_file`

---

## Result Structure

All functions return a dictionary with a `status` field:
- `"success"`: Operation completed successfully
- `"error"`: Operation failed (check `message` field for details)

## Example

- [multi_layer_filling.R](multi_layer_filling.R) - Multi-layer CO2 filling with synthetic geology (R equivalent of `examples/multi_layer_filling.jl`)

## Notes

- **Grid indexing**: R uses 1-based indexing, which matches the Julia interface
- **Layer ordering**: Pass layers shallowest-first, deepest-last. The simulator reverses them internally.
- **Boundary padding**: With closed boundaries and `pad_width=1` (default), injection matrices must have dimensions `(n_times, nx + 2, ny + 2)`. The `nx_after_bc` and `ny_after_bc` values from `setup_simulator_from_surfaces` give the correct padded size.
- **Array dimensions**: Injection matrices use dimension order `(n_times x nx_after_bc x ny_after_bc)`
- **Units**:
  - Injection rates: m^3/year
  - Volumes: m^3
  - Time: years
  - Pressure: Pa
  - Density: kg/m^3
- **Performance**: First simulation run will be slower due to Julia JIT compilation. Subsequent runs are much faster.
