"""
Monte Carlo simulation runner for CO2 injection uncertainty quantification.

This module handles:
- Generating parameter samples from distributions
- Running simulations with sampled parameters
- Collecting and aggregating results
"""

using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using Random
using Distributions
using ProgressMeter

include("monte_carlo_config.jl")

"""
    MonteCarloResult

Stores results from a single Monte Carlo realization.
"""
struct MonteCarloResult
    realization_id::Int
    parameters::Dict{Symbol, Float64}
    snapshots::Vector{ReservoirSnapshot}
end

"""
    MonteCarloEnsemble

Aggregated results from all Monte Carlo realizations.
"""
struct MonteCarloEnsemble
    config::MonteCarloConfig
    results::Vector{MonteCarloResult}
    timepoints::Vector{Float64}
end

"""
    sample_reservoir_properties(config::MonteCarloConfig, n_layers::Int)

Generate a random sample of reservoir properties for all layers.

# Arguments
- `config::MonteCarloConfig`: Monte Carlo configuration with parameter distributions
- `n_layers::Int`: Number of reservoir layers

# Returns
- `Vector{ReservoirProperties}`: Sampled properties for each layer
"""
function sample_reservoir_properties(config::MonteCarloConfig, n_layers::Int)
    # Sample parameters from distributions
    sampled_params = Dict{Symbol, Float64}()
    for (param_name, dist) in config.param_distributions
        sampled_params[param_name] = rand(dist)
    end

    # Density values from L1 up to L9 (from paper)
    brine_density = 1020.0
    co2_density = (570.0 + 350.0) / 2.0  # Average density between L1 and L9

    # Get spatial variability parameter for pressure threshold (default to 0.0 if not specified)
    shale_pressure_threshold_std = get(sampled_params, :shale_pressure_threshold_std, 0.0)

    # Create ReservoirProperties for each layer
    # Note: leakage_height is computed automatically from shale_pressure_threshold
    reservoir_properties = Vector{ReservoirProperties}(undef, n_layers)

    for i in 1:n_layers
        # Top layer has infinite leakage height (impermeable caprock)
        if i == n_layers
            reservoir_properties[i] = ReservoirProperties(
                sampled_params[:sand_porosity],
                sampled_params[:sand_residual_co2_saturation],
                sampled_params[:sand_irreducible_water_saturation],
                sampled_params[:shale_pressure_threshold],
                sampled_params[:residual_leakage_time];
                leakage_height=Inf,  # Explicitly set to Inf for top layer
                shale_pressure_threshold_std=shale_pressure_threshold_std,
                brine_density=brine_density,
                co2_density=co2_density
            )
        else
            reservoir_properties[i] = ReservoirProperties(
                sampled_params[:sand_porosity],
                sampled_params[:sand_residual_co2_saturation],
                sampled_params[:sand_irreducible_water_saturation],
                sampled_params[:shale_pressure_threshold],
                sampled_params[:residual_leakage_time];
                # leakage_height computed automatically from pressure
                shale_pressure_threshold_std=shale_pressure_threshold_std,
                brine_density=brine_density,
                co2_density=co2_density
            )
        end
    end

    return reservoir_properties, sampled_params
end

"""
    run_single_realization(
        realization_id::Int,
        config::MonteCarloConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        injection_events::Vector{Vector{InjectionEvent}}
    )

Run a single Monte Carlo realization with sampled parameters.

# Arguments
- `realization_id::Int`: Identifier for this realization
- `config::MonteCarloConfig`: Monte Carlo configuration
- `layers::Vector{Layer}`: Reservoir layer structure
- `domain::Domain3D`: Spatial domain
- `injection_events::Vector{Vector{InjectionEvent}}`: CO2 injection schedule (per layer)

# Returns
- `MonteCarloResult`: Results from this realization
"""
function run_single_realization(
    realization_id::Int,
    config::MonteCarloConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    injection_events::Vector{Vector{InjectionEvent}}
)
    # Sample reservoir properties
    n_layers = length(layers)
    reservoir_properties, sampled_params = sample_reservoir_properties(config, n_layers)

    # Run simulation
    seqs, leakage_states = fill_layers(
        layers,
        domain,
        reservoir_properties,
        injection_events;
        verbose=false
    )

    # Generate snapshots
    snapshots = generate_reservoir_snapshots(
        layers,
        seqs,
        leakage_states,
        domain,
        reservoir_properties,
        injection_events;
        num_snapshots=config.num_snapshots,
        start_time=config.start_time,
        end_time=config.end_time,
        verbose=false
    )

    return MonteCarloResult(realization_id, sampled_params, snapshots)
end

"""
    run_monte_carlo_simulation(
        config::MonteCarloConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        injection_events::Vector{Vector{InjectionEvent}};
        verbose::Bool=true
    )

Run Monte Carlo uncertainty analysis with multiple parameter realizations.

# Arguments
- `config::MonteCarloConfig`: Monte Carlo configuration
- `layers::Vector{Layer}`: Reservoir layer structure
- `domain::Domain3D`: Spatial domain
- `injection_events::Vector{Vector{InjectionEvent}}`: CO2 injection schedule (per layer)
- `verbose::Bool=true`: Show progress bar

# Returns
- `MonteCarloEnsemble`: Aggregated results from all realizations
"""
function run_monte_carlo_simulation(
    config::MonteCarloConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    injection_events::Vector{Vector{InjectionEvent}};
    verbose::Bool=true
)
    # Set random seed if specified
    if config.random_seed !== nothing
        Random.seed!(config.random_seed)
    end

    # Compute timepoints
    timepoints = range(config.start_time, config.end_time, length=config.num_snapshots)

    # Initialize results storage
    results = Vector{MonteCarloResult}(undef, config.n_realizations)

    # Progress bar
    progress = Progress(config.n_realizations, desc="Running Monte Carlo simulations: ", enabled=verbose)

    # Run simulations
    for i in 1:config.n_realizations
        results[i] = run_single_realization(i, config, layers, domain, injection_events)
        next!(progress)
    end

    return MonteCarloEnsemble(config, results, collect(timepoints))
end

"""
    extract_timeseries(ensemble::MonteCarloEnsemble, field::Symbol)

Extract a specific field from all realizations as a matrix.

# Arguments
- `ensemble::MonteCarloEnsemble`: Monte Carlo ensemble results
- `field::Symbol`: Field to extract (e.g., :total_stored_m3, :total_leaked_m3)

# Returns
- `Matrix{Float64}`: Matrix of size (n_timepoints, n_realizations)
"""
function extract_timeseries(ensemble::MonteCarloEnsemble, field::Symbol)
    n_timepoints = length(ensemble.timepoints)
    n_realizations = length(ensemble.results)

    data = Matrix{Float64}(undef, n_timepoints, n_realizations)

    for (i, result) in enumerate(ensemble.results)
        for (j, snapshot) in enumerate(result.snapshots)
            data[j, i] = getfield(snapshot, field)
        end
    end

    return data
end

"""
    extract_layer_timeseries(ensemble::MonteCarloEnsemble, layer_idx::Int)

Extract CO2 volume in a specific layer for all realizations.

# Arguments
- `ensemble::MonteCarloEnsemble`: Monte Carlo ensemble results
- `layer_idx::Int`: Layer index (1-based)

# Returns
- `Matrix{Float64}`: Matrix of size (n_timepoints, n_realizations)
"""
function extract_layer_timeseries(ensemble::MonteCarloEnsemble, layer_idx::Int)
    n_timepoints = length(ensemble.timepoints)
    n_realizations = length(ensemble.results)

    data = Matrix{Float64}(undef, n_timepoints, n_realizations)

    for (i, result) in enumerate(ensemble.results)
        for (j, snapshot) in enumerate(result.snapshots)
            data[j, i] = snapshot.stored_by_layer_m3[layer_idx]
        end
    end

    return data
end

"""
    compute_percentiles(data::Matrix{Float64}, percentiles::Vector{Float64})

Compute percentiles along realizations (columns) for each timepoint (row).

# Arguments
- `data::Matrix{Float64}`: Data matrix (n_timepoints × n_realizations)
- `percentiles::Vector{Float64}`: Percentiles to compute (e.g., [5, 50, 95])

# Returns
- `Matrix{Float64}`: Percentile matrix (n_timepoints × n_percentiles)
"""
function compute_percentiles(data::Matrix{Float64}, percentiles::Vector{Float64})
    n_timepoints = size(data, 1)
    n_percentiles = length(percentiles)

    result = Matrix{Float64}(undef, n_timepoints, n_percentiles)

    for i in 1:n_timepoints
        for (j, p) in enumerate(percentiles)
            result[i, j] = quantile(data[i, :], p / 100.0)
        end
    end

    return result
end
