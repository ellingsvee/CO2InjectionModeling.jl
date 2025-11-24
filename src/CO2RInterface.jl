module CO2RInterface

using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using Interpolations

export setup_simulator, configure_reservoir, run_simulation

# Global state to hold simulator configuration
mutable struct SimulatorState
    topography::Union{Nothing, SleipnerTopography}
    domain::Union{Nothing, Domain3D}
    layers::Union{Nothing, Vector{Layer}}
    reservoir_properties::Union{Nothing, Vector{ReservoirProperties}}
    boundary_condition::Symbol
end

const SIMULATOR = SimulatorState(nothing, nothing, nothing, nothing, :open)

"""
    setup_simulator(; data_path="sleipner/depth_surfaces/", boundary_condition="open")

Load and set up the Sleipner topography layers for simulation.

# Arguments
- `data_path`: Path to the Sleipner depth surfaces data (default: "sleipner/depth_surfaces/")
- `boundary_condition`: Boundary condition type, either "open" or "closed" (default: "open")

# Returns
- Dictionary with setup information:
  - `n_layers`: Number of layers
  - `nx`, `ny`: Grid dimensions for each layer
  - `status`: "success" or error message

# Example (R)
```
library(JuliaCall)
julia_setup()
julia_source("path/to/CO2InjectionModeling.jl")

# Setup the simulator
result <- julia_call("setup_simulator",
                     data_path = "sleipner/depth_surfaces/",
                     boundary_condition = "open")
print(result)
```
"""
function setup_simulator(;
    data_path::String="sleipner/depth_surfaces/",
    boundary_condition::String="open"
)
    try
        # Validate boundary condition
        bc_symbol = Symbol(boundary_condition)
        if !(bc_symbol in [:open, :closed])
            return Dict(
                "status" => "error",
                "message" => "boundary_condition must be 'open' or 'closed'"
            )
        end

        # Load topography
        SIMULATOR.topography = load_sleipner_topography(data_path)
        SIMULATOR.domain = create_domain_from_topography(SIMULATOR.topography, 1.0)
        SIMULATOR.layers = analyze_base_surfaces(SIMULATOR.topography; boundary_condition=bc_symbol)
        SIMULATOR.boundary_condition = bc_symbol

        # Initialize with default reservoir properties
        SIMULATOR.reservoir_properties = generate_reservoir_properties_for_sleipner_layers()

        return Dict(
            "status" => "success",
            "n_layers" => length(SIMULATOR.layers),
            "nx" => SIMULATOR.topography.nx,
            "ny" => SIMULATOR.topography.ny,
            "boundary_condition" => boundary_condition
        )
    catch e
        return Dict(
            "status" => "error",
            "message" => string(e)
        )
    end
end

"""
    configure_reservoir(porosity, residual_co2_sat, irreducible_water_sat,
                       shale_pressure_threshold, brine_co2_density_diff,
                       residual_leakage_time; layer_specific=false)

Configure reservoir properties for the simulation.

# Arguments
- `porosity`: Sand porosity (0-1). Can be single value or vector of length n_layers
- `residual_co2_sat`: Residual CO2 saturation (0-1). Can be single value or vector
- `irreducible_water_sat`: Irreducible water saturation (0-1). Can be single value or vector
- `shale_pressure_threshold`: Shale pressure threshold (Pa). Can be single value or vector
- `brine_co2_density_diff`: Density difference between brine and CO2 (kg/m³). Can be single value or vector
- `residual_leakage_time`: Residual leakage time (years). Can be single value or vector
- `layer_specific`: If true, expects vectors for layer-specific properties (default: false)

# Returns
- Dictionary with configuration status

# Example (R)
```
# Use same properties for all layers
result <- julia_call("configure_reservoir",
                     porosity = 0.35,
                     residual_co2_sat = 0.2,
                     irreducible_water_sat = 0.3,
                     shale_pressure_threshold = 98000.0,
                     brine_co2_density_diff = 450.0,
                     residual_leakage_time = 1.0,
                     layer_specific = FALSE)

# Or specify different properties for each layer (9 layers for Sleipner)
result <- julia_call("configure_reservoir",
                     porosity = c(0.35, 0.36, 0.37, 0.35, 0.36, 0.37, 0.35, 0.36, 0.37),
                     residual_co2_sat = rep(0.2, 9),
                     irreducible_water_sat = rep(0.3, 9),
                     shale_pressure_threshold = rep(98000.0, 9),
                     brine_co2_density_diff = c(450, 477.5, 505, 532.5, 560, 587.5, 615, 642.5, 670),
                     residual_leakage_time = rep(1.0, 9),
                     layer_specific = TRUE)
```
"""
function configure_reservoir(
    porosity::Union{Float64, Vector{Float64}},
    residual_co2_sat::Union{Float64, Vector{Float64}},
    irreducible_water_sat::Union{Float64, Vector{Float64}},
    shale_pressure_threshold::Union{Float64, Vector{Float64}},
    brine_co2_density_diff::Union{Float64, Vector{Float64}},
    residual_leakage_time::Union{Float64, Vector{Float64}};
    layer_specific::Bool=false
)
    try
        if isnothing(SIMULATOR.layers)
            return Dict(
                "status" => "error",
                "message" => "Must call setup_simulator first"
            )
        end

        n_layers = length(SIMULATOR.layers)

        # Convert scalars to vectors if not layer_specific
        if !layer_specific
            porosity = fill(Float64(porosity), n_layers)
            residual_co2_sat = fill(Float64(residual_co2_sat), n_layers)
            irreducible_water_sat = fill(Float64(irreducible_water_sat), n_layers)
            shale_pressure_threshold = fill(Float64(shale_pressure_threshold), n_layers)
            brine_co2_density_diff = fill(Float64(brine_co2_density_diff), n_layers)
            residual_leakage_time = fill(Float64(residual_leakage_time), n_layers)
        end

        # Validate vector lengths
        if length(porosity) != n_layers ||
           length(residual_co2_sat) != n_layers ||
           length(irreducible_water_sat) != n_layers ||
           length(shale_pressure_threshold) != n_layers ||
           length(brine_co2_density_diff) != n_layers ||
           length(residual_leakage_time) != n_layers
            return Dict(
                "status" => "error",
                "message" => "All property vectors must have length $n_layers"
            )
        end

        # Create ReservoirProperties for each layer
        SIMULATOR.reservoir_properties = [
            ReservoirProperties(
                porosity[i],
                residual_co2_sat[i],
                irreducible_water_sat[i],
                shale_pressure_threshold[i],
                brine_co2_density_diff[i],
                residual_leakage_time[i]
            )
            for i in 1:n_layers
        ]

        return Dict(
            "status" => "success",
            "n_layers" => n_layers
        )
    catch e
        return Dict(
            "status" => "error",
            "message" => string(e)
        )
    end
end

"""
    run_simulation(start_time, end_time, time_step, injection_times,
                  injection_locations_i, injection_locations_j,
                  injection_amounts, injection_layer_indices;
                  num_snapshots=nothing, verbose=false)

Run a CO2 injection simulation with specified parameters.

# Arguments
- `start_time`: Simulation start time (years)
- `end_time`: Simulation end time (years)
- `time_step`: Time step between injection events (years)
- `injection_times`: Vector of injection event times (years)
- `injection_locations_i`: Vector of i-indices for injection locations (1-based)
- `injection_locations_j`: Vector of j-indices for injection locations (1-based)
- `injection_amounts`: Vector of injection amounts (m³/year at each location/time)
- `injection_layer_indices`: Vector of layer indices for each injection (1-based, 1=bottom)
- `num_snapshots`: Number of snapshots to save (default: auto-calculated from time_step)
- `verbose`: Print progress messages (default: false)

# Returns
- Dictionary containing simulation summary with:
  - `timepoints`: Vector of snapshot times
  - `total_co2_volumes`: Total CO2 volume at each timepoint
  - `layer_co2_volumes`: Matrix of CO2 volumes per layer (timepoints × layers)
  - `trap_co2_volumes`: List of matrices for trap volumes in each layer
  - `trap_co2_percentages`: List of matrices for trap percentages in each layer
  - `num_layers`: Number of layers
  - `num_traps_per_layer`: Vector of trap counts per layer

# Example (R)
```
# Simple example: inject at one location in bottom layer
result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 15.0,
                     time_step = 1.0,
                     injection_times = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                     injection_locations_i = rep(32L, 15),
                     injection_locations_j = rep(59L, 15),
                     injection_amounts = c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
                                          0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
                                          0.82, 0.86, 0.76) * 1e9 / 570,
                     injection_layer_indices = rep(1L, 15),
                     num_snapshots = 14L,
                     verbose = TRUE)

# Access results
print(result\$timepoints)
print(result\$total_co2_volumes)
```
"""
function run_simulation(
    start_time::Float64,
    end_time::Float64,
    time_step::Float64,
    injection_times::Vector{Float64},
    injection_locations_i::Vector{Int},
    injection_locations_j::Vector{Int},
    injection_amounts::Vector{Float64},
    injection_layer_indices::Vector{Int};
    num_snapshots::Union{Nothing, Int}=nothing,
    verbose::Bool=false
)
    try
        if isnothing(SIMULATOR.layers) || isnothing(SIMULATOR.reservoir_properties)
            return Dict(
                "status" => "error",
                "message" => "Must call setup_simulator and configure_reservoir first"
            )
        end

        # Validate inputs
        n_injections = length(injection_times)
        if length(injection_locations_i) != n_injections ||
           length(injection_locations_j) != n_injections ||
           length(injection_amounts) != n_injections ||
           length(injection_layer_indices) != n_injections
            return Dict(
                "status" => "error",
                "message" => "All injection vectors must have the same length"
            )
        end

        # Calculate num_snapshots if not provided
        if isnothing(num_snapshots)
            num_snapshots = Int(ceil((end_time - start_time) / time_step))
        end

        # Create injection events for each layer
        n_layers = length(SIMULATOR.layers)
        grid_size = size(SIMULATOR.layers[1].trap_structure.topography)

        # Group injections by layer and time
        injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)

        for layer_idx in 1:n_layers
            # Find all injections for this layer
            layer_mask = injection_layer_indices .== layer_idx
            layer_times = injection_times[layer_mask]
            layer_i = injection_locations_i[layer_mask]
            layer_j = injection_locations_j[layer_mask]
            layer_amounts = injection_amounts[layer_mask]

            if isempty(layer_times)
                # No injection in this layer
                injection_events[layer_idx] = [InjectionEvent(0.0, zeros(grid_size))]
            else
                # Create injection events for this layer
                unique_times = unique(layer_times)
                sort!(unique_times)

                layer_events = InjectionEvent[]
                for t in unique_times
                    time_mask = layer_times .== t
                    injection_rate = zeros(grid_size)

                    # Add all injections at this time
                    for idx in findall(time_mask)
                        i = layer_i[idx]
                        j = layer_j[idx]
                        amount = layer_amounts[idx]
                        injection_rate[i, j] += amount
                    end

                    push!(layer_events, InjectionEvent(t, injection_rate))
                end

                injection_events[layer_idx] = layer_events
            end
        end

        # Run simulation
        if verbose
            println("Running simulation from $start_time to $end_time years with $num_snapshots snapshots")
        end

        seqs, leakage_states, snapshots = fill_layers(
            SIMULATOR.layers,
            SIMULATOR.domain,
            SIMULATOR.reservoir_properties,
            injection_events;
            num_snapshots=num_snapshots,
            start_time=start_time,
            end_time=end_time,
            verbose=verbose
        )

        if verbose
            println("Generating simulation summary...")
        end

        summary = generate_simulation_summary(snapshots, SIMULATOR.layers, seqs)

        # Convert summary to Dict for R
        result = Dict(
            "status" => "success",
            "timepoints" => summary.timepoints,
            "total_co2_volumes" => summary.total_co2_volumes,
            "layer_co2_volumes" => summary.layer_co2_volumes,
            "trap_co2_volumes" => summary.trap_co2_volumes,
            "trap_co2_percentages" => summary.trap_co2_percentages,
            "num_layers" => summary.num_layers,
            "num_traps_per_layer" => summary.num_traps_per_layer
        )

        if verbose
            println("Simulation completed successfully")
        end

        return result
    catch e
        return Dict(
            "status" => "error",
            "message" => string(e),
            "stacktrace" => sprint(showerror, e, catch_backtrace())
        )
    end
end

end # module
