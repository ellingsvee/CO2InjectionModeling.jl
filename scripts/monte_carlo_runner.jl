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
    seqs::Union{Vector{Vector{SpillEvent}}, Nothing}
    leakage_states::Union{Vector{LeakageState}, Nothing}
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

    # Special handling for the shale_pressure_threshold as we need to sample per layer
    sampled_shale_thresholds = rand.(Ref(config.param_distributions[:shale_pressure_threshold]), n_layers)
    for i in 1:n_layers
        sampled_params[Symbol("shale_pressure_threshold_layer_$i")] = sampled_shale_thresholds[i]
    end

    # Density values from L1 up to L9 (from paper)
    brine_density = 1020.0
    co2_density = (570.0 + 350.0) / 2.0  # Average density between L1 and L9

    # Get spatial variability parameter for pressure threshold (default to 0.0 if not specified)
    shale_pressure_threshold_std = get(sampled_params, :shale_pressure_threshold_std, 0.0)

    # Create ReservoirProperties for each layer
    # Note: leakage_height is computed automatically from shale_pressure_threshold
    reservoir_properties = Vector{ReservoirProperties}(undef, n_layers)


    # sampled_shale_thresholds =[
    #     64935.189136968474,
    #     145480.99151323555,
    #     60265.54535253754,
    #     86545.36470224884,
    #     134536.57733389136,
    #     114664.10814779482,
    #     118790.15494120165,
    #     110065.36326170494,
    #     115946.26594344187,
    #     # Inf, 
    # ]
    # sampled_shale_thresholds[1] *= 1.5
    # sampled_shale_thresholds[2] *= 0.95
    # sampled_shale_thresholds[3] *= 1.75
    # sampled_shale_thresholds[4] *= 1.25

    # sampled_shale_thresholds[6] *= 0.88
    # sampled_shale_thresholds[7] *= 0.85
    # sampled_shale_thresholds[8] *= 1.1
    # sampled_shale_thresholds[9] *= 1.02

    sampled_shale_thresholds =[
        64935.189136968474,
        145480.99151323555,
        60265.54535253754,
        86545.36470224884,
        134536.57733389136,
        114664.10814779482,
        118790.15494120165,
        110065.36326170494,
        115946.26594344187,
        # Inf, 
    ]
    sampled_shale_thresholds[1] *= 1.5
    sampled_shale_thresholds[2] *= 0.6
    sampled_shale_thresholds[3] *= 1.75
    sampled_shale_thresholds[4] *= 1.25
    sampled_shale_thresholds[5] *= 1.1

    sampled_shale_thresholds[6] *= 0.88
    sampled_shale_thresholds[7] *= 0.85
    sampled_shale_thresholds[8] *= 1.1
    sampled_shale_thresholds[9] *= 1.0

    for i in 1:n_layers
        # Top layer has infinite leakage height (impermeable caprock)
        # if i == n_layers
        #     sampled_params[Symbol("shale_pressure_threshold_layer_$i")] = Inf
        # end

        reservoir_properties[i] = ReservoirProperties(
            sampled_params[:sand_porosity],
            sampled_params[:sand_residual_co2_saturation],
            sampled_params[:sand_irreducible_water_saturation],
            # sampled_params[Symbol("shale_pressure_threshold_layer_$i")],
            sampled_shale_thresholds[i],
            sampled_params[:residual_leakage_time];
            # leakage_height computed automatically from pressure
            shale_pressure_threshold_std=shale_pressure_threshold_std,
            brine_density=brine_density,
            co2_density=co2_density
        )

    end

    return reservoir_properties, sampled_params
end

"""
    run_single_realization(
        realization_id::Int,
        config::MonteCarloConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        injection_events::Vector{Vector{InjectionEvent}};
        store_seqs::Bool=false
    )

Run a single Monte Carlo realization with sampled parameters.

# Arguments
- `realization_id::Int`: Identifier for this realization
- `config::MonteCarloConfig`: Monte Carlo configuration
- `layers::Vector{Layer}`: Reservoir layer structure
- `domain::Domain3D`: Spatial domain
- `injection_events::Vector{Vector{InjectionEvent}}`: CO2 injection schedule (per layer)
- `store_seqs::Bool=false`: Whether to store spill sequences for spatial animations

# Returns
- `MonteCarloResult`: Results from this realization
"""
function run_single_realization(
    realization_id::Int,
    config::MonteCarloConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    injection_events::Vector{Vector{InjectionEvent}};
    store_seqs::Bool=false
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

    # Optionally store seqs and leakage_states for spatial animations
    stored_seqs = store_seqs ? seqs : nothing
    stored_leakage_states = store_seqs ? leakage_states : nothing

    return MonteCarloResult(realization_id, sampled_params, snapshots, stored_seqs, stored_leakage_states)
end

"""
    run_monte_carlo_simulation(
        config::MonteCarloConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        injection_events::Vector{Vector{InjectionEvent}};
        verbose::Bool=true,
        store_seqs::Bool=false
    )

Run Monte Carlo uncertainty analysis with multiple parameter realizations.

# Arguments
- `config::MonteCarloConfig`: Monte Carlo configuration
- `layers::Vector{Layer}`: Reservoir layer structure
- `domain::Domain3D`: Spatial domain
- `injection_events::Vector{Vector{InjectionEvent}}`: CO2 injection schedule (per layer)
- `verbose::Bool=true`: Show progress bar
- `store_seqs::Bool=false`: Whether to store spill sequences for spatial animations

# Returns
- `MonteCarloEnsemble`: Aggregated results from all realizations
"""
function run_monte_carlo_simulation(
    config::MonteCarloConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    injection_events::Vector{Vector{InjectionEvent}};
    verbose::Bool=true,
    store_seqs::Bool=false
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
        results[i] = run_single_realization(i, config, layers, domain, injection_events; store_seqs=store_seqs)
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
