# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

CO2InjectionModeling.jl simulates CO2 injection in multi-layered subsurface geological formations using invasion percolation modeling. The primary application is modeling the Sleipner CO2 storage site with 9 sandstone layers separated by shale barriers. The code uses the SurfaceWaterIntegratedModeling (SWIM) library for trap-based flow simulation.

## Common Commands

### Setup and Installation

```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Activate the project environment in Julia REPL
julia --project=.
```

### Running Simulations

```bash
# Run single-layer test
julia --project=. scripts/test_single_layer.jl

# Run multi-layer test
julia --project=. scripts/test_multi_layer.jl

# Run Monte Carlo uncertainty analysis
julia --project=. scripts/monte_carlo_simulation.jl

# Find optimal parameters for Sleipner match
julia --project=. scripts/find_optimal_params.jl
```

### Testing

There is currently no formal test suite. Scripts in `scripts/` serve as integration tests and examples.

## Architecture

### Core Module Structure

The main module (`src/CO2InjectionModeling.jl`) includes submodules in a specific order due to dependencies:

1. **sleipner_setup.jl** - Loads Sleipner topography data from .npy files, defines injection schedules, and converts UTM coordinates to grid indices
2. **structs.jl** - Core data structures (Domain3D, ReservoirProperties, LeakageState, etc.)
3. **layer_analysis.jl** - Analyzes trap structures for each geological layer using SWIM
4. **volume_conversion.jl** - Converts between SWIM units and physical units (m³)
5. **utils.jl** - Utilities including trap hierarchy traversal (must load before leakage.jl)
6. **leakage.jl** - Leakage detection, residual drainage, and inter-layer flow modeling
7. **fill_layer.jl** - Single-layer filling simulation
8. **fill_layers.jl** - Multi-layer simulation orchestration with layer-to-layer leakage
9. **analysis.jl** - Post-processing: generates snapshots, computes mass balance, statistics
10. **visualization.jl** - Animation generation (bird's eye view of CO2 distribution)
11. **CO2RInterface.jl** - R interface via JuliaCall for external access

### Key Data Flow

1. **Topography Loading**: `load_sleipner_topography()` reads .npy depth surfaces from `sleipner/depth_surfaces/`, reverses Y-axis for UTM convention
2. **Layer Analysis**: `analyze_base_surfaces()` creates `TrapStructure` for each layer using SWIM's spill analysis
3. **Injection Events**: `generate_sleipner_injection_events()` creates time-varying injection schedule (1996-2010 historical data)
4. **Simulation**: `fill_layers()` orchestrates multi-layer filling:
   - Bottom layer receives injection from well
   - When CO2 column height exceeds `leakage_height` in a trap, CO2 migrates upward
   - Residual drainage modeled with `residual_leakage_time`
   - Each layer's output becomes `WeatherEvent` input for the layer above
5. **Analysis**: `generate_reservoir_snapshots()` post-processes spill sequences into physical metrics (volumes, heights, mass balance)

### Leakage Modeling

Leakage between layers is complex and performance-critical:

- **Pressure-based thresholds**: `shale_pressure_threshold` (Pa) converted to `leakage_height` (m) via buoyancy: `h = P / (Δρ * g)`
- **Per-trap sampling**: Each trap can have different leakage height (spatial heterogeneity via `shale_pressure_threshold_std`)
- **Residual drainage**: After leakage starts, trapped CO2 drains over time with residual saturation (`sand_residual_co2_saturation`)
- **Leak location**: Identified as the lowest spill point in the trap (see `find_leakage_location()`)
- **Ancestor drainage**: When a trap leaks, CO2 in filled parent traps also drains through it
- **WeatherEvent generation**: Leakage from layer N becomes time-varying input to layer N+1

### Monte Carlo Uncertainty Analysis

Located in `scripts/monte_carlo_*.jl`:

- **Config** (`monte_carlo_config.jl`): Defines parameter distributions (Normal, Uniform, etc.) for uncertain parameters
- **Runner** (`monte_carlo_runner.jl`): Samples parameters, runs multiple realizations, collects results
- **Simulation** (`monte_carlo_simulation.jl`): High-level script that sets up and runs ensemble
- **Analysis** (`monte_carlo_analysis.jl`): Computes percentiles, generates uncertainty plots
- **Parameter fitting** (`find_optimal_params.jl`): Optimizes shale pressure thresholds to match seismic observations

Typical workflow: Sample `shale_pressure_threshold` per layer, run 100-500 realizations, compare mass distribution across layers to seismic data.

## Sleipner-Specific Information

### Grid and Layers

- **Grid**: 64 × 118 horizontal cells, ~50m × 50m resolution, 3.2 km × 5.9 km extent
- **9 Sand Layers**: L1 (deepest, injection layer) to L9 (shallowest, below caprock)
- **8 Shale Barriers**: 7 thin intra-formational shales (~1.5m) + 1 thick shale between L8/L9 (~7.5m)
- **Injection Point**: Grid cell approximately (32, 59), located in L1, 65m north of main feeder chimney
- **Topography Orientation**: Y-axis is reversed in `load_sleipner_topography()` to match UTM northing convention

### Physical Parameters (from `sleipner_physical_quanities.md`)

- **Porosity**: 0.4 (sandstone)
- **Residual CO2 saturation**: 0.2
- **Irreducible water saturation**: 0.3
- **Shale pressure threshold**: 98 kPa (base case), varied 49-147 kPa in sensitivity analysis
- **CO2 density**: Depth-dependent, 350 kg/m³ (L9) to 570 kg/m³ (L1); code currently uses average ~460 kg/m³
- **Brine density**: 1020 kg/m³
- **Injection period**: 1996-2010, total 12.18 Mt (see `generate_sleipner_injection_events()` for annual rates)

### Feeder Chimneys

Vertical conduits through shale layers (relevant for "Shales with Breaks" scenarios):
- **Main feeder**: ~100m × 100m, 65m south of injection point, parsed from Z-MAP polygon in `sleipner/feeders/data/Main_feeder_chimney`
- Use `load_feeder_location()` to get centroid and grid index
- Use `utm_to_grid_index()` to convert arbitrary UTM coordinates to grid cells

## Important Implementation Notes

### Coordinate Systems

- **UTM coordinates**: Real-world easting/northing (meters), origin ~(436800, 6468100)
- **Grid indices**: 1-based Julia indexing (i, j), origin at southwest corner
- **Y-axis convention**: Array indices increase downward, but UTM northing increases upward → topography surfaces are reversed on load

### Performance Considerations

- Leakage detection runs every time step for all traps → use cached interpolations (`CachedInterpolations`) to avoid repeated volume-to-height conversions
- Monte Carlo simulations are embarrassingly parallel but currently run serially → consider parallelization for large ensembles
- SWIM's trap analysis is one-time setup cost, but fill sequence computation scales with number of injection events

### Common Pitfalls

- **Top layer (L9)**: Must have `leakage_height = Inf` to represent impermeable caprock
- **Layer ordering**: L1 is deepest (index 1), L9 is shallowest (index 9) in code arrays
- **Volume units**: SWIM uses internal units; always convert to m³ using functions in `volume_conversion.jl`
- **Injection matrix dimensions**: When using R interface, injection arrays are (n_times × nx × ny)
- **Boundary conditions**: "open" allows lateral spillage; "closed" adds 1-cell border (affects effective grid size)

## R Interface

The package can be called from R via JuliaCall (see `r_interface/R_INTERFACE_README.md`):

1. `setup_simulator()` - Load topography and initialize
2. `setup_sleipner_reservoir()` or `configure_reservoir()` - Set reservoir properties
3. `run_simulation()` - Execute simulation with injection schedule
4. `generate_birdseye_animation()` - Create visualization

Key functions are exported from `CO2RInterface` module with global state in `SIMULATOR`.

## Code Style

- Functions use lowercase with underscores (e.g., `load_sleipner_topography`)
- Types use CamelCase (e.g., `ReservoirProperties`)
- Exported symbols are explicitly listed at the top of each file
- Verbose mode controlled by `verbose::Bool` parameter (default false)
- Use `@warn` or `println` for user-facing messages; avoid excessive debug output
