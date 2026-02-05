# Injection schedule optimization for multi-well CO2 injection
# Maximizes CO2 storage by optimizing time-varying injection rates at multiple wells

using Random
using BlackBoxOptim

export OptimizationConfig, BayesianOptimizationConfig, InjectionOptimizationResult
export optimize_injection_schedule, evaluate_injection_schedule
export create_multi_well_injection_events
export optimize_sleipner_injection
export plot_injection_schedule, plot_optimization_convergence

# Location optimization exports
export LocationOptimizationConfig, LocationOptimizationResult
export optimize_well_locations, evaluate_well_locations
export optimize_sleipner_well_locations
export plot_well_locations

# =============================================================================
# DATA STRUCTURES
# =============================================================================

"""
    OptimizationConfig

Configuration for injection schedule optimization.

# Fields
- `n_time_periods::Int`: Number of time periods to divide the injection horizon
- `well_locations::Vector{CartesianIndex}`: Grid indices for each injection well
- `well_layers::Vector{Int}`: Layer index for each well (1 = bottom/L1, 9 = top/L9)
- `total_mass_mt::Float64`: Total CO2 mass to inject (Mt)
- `max_rate_mt_per_year::Float64`: Maximum injection rate per well (Mt/year)
- `min_rate_mt_per_year::Float64`: Minimum injection rate per well (Mt/year)
- `start_time::Float64`: Simulation start time (years)
- `end_time::Float64`: Simulation end time (years)
- `co2_density::Float64`: CO2 density for volume conversion (kg/m³)
"""
struct OptimizationConfig
    n_time_periods::Int
    well_locations::Vector{CartesianIndex{2}}
    well_layers::Vector{Int}
    total_mass_mt::Float64
    max_rate_mt_per_year::Float64
    min_rate_mt_per_year::Float64
    start_time::Float64
    end_time::Float64
    co2_density::Float64

    function OptimizationConfig(
        n_time_periods::Int,
        well_locations::Vector{<:CartesianIndex},
        well_layers::Vector{Int},
        total_mass_mt::Float64,
        max_rate_mt_per_year::Float64,
        min_rate_mt_per_year::Float64,
        start_time::Float64,
        end_time::Float64,
        co2_density::Float64
    )
        @assert n_time_periods >= 1 "Need at least 1 time period"
        @assert length(well_locations) == length(well_layers) "Well locations and layers must match"
        @assert length(well_locations) >= 1 "Need at least 1 well"
        @assert max_rate_mt_per_year >= min_rate_mt_per_year "Max rate must be >= min rate"
        @assert end_time > start_time "End time must be after start time"
        @assert total_mass_mt > 0 "Total mass must be positive"

        new(n_time_periods, collect(CartesianIndex{2}, well_locations), well_layers, total_mass_mt,
            max_rate_mt_per_year, min_rate_mt_per_year, start_time, end_time, co2_density)
    end
end

# Convenience constructor when all wells are in layer 1
function OptimizationConfig(
    n_time_periods::Int,
    well_locations::Vector{<:CartesianIndex},
    total_mass_mt::Float64,
    max_rate_mt_per_year::Float64,
    min_rate_mt_per_year::Float64,
    start_time::Float64,
    end_time::Float64,
    co2_density::Float64
)
    well_layers = fill(1, length(well_locations))
    return OptimizationConfig(
        n_time_periods, well_locations, well_layers, total_mass_mt,
        max_rate_mt_per_year, min_rate_mt_per_year, start_time, end_time, co2_density
    )
end


"""
    BayesianOptimizationConfig

Configuration for Bayesian Optimization algorithm.

# Fields
- `n_initial_samples::Int`: Number of initial random samples for surrogate training
- `max_evaluations::Int`: Maximum objective function evaluations
- `convergence_threshold::Float64`: Stop if improvement < threshold for N iterations
- `convergence_patience::Int`: Number of iterations to wait before early stopping
"""
struct BayesianOptimizationConfig
    n_initial_samples::Int
    max_evaluations::Int
    convergence_threshold::Float64
    convergence_patience::Int

    function BayesianOptimizationConfig(;
        n_initial_samples::Int = 20,
        max_evaluations::Int = 100,
        convergence_threshold::Float64 = 1e-4,
        convergence_patience::Int = 10
    )
        new(n_initial_samples, max_evaluations, convergence_threshold, convergence_patience)
    end
end


"""
    InjectionOptimizationResult

Result from injection schedule optimization.

# Fields
- `optimal_rates::Matrix{Float64}`: Optimal injection rates (n_periods x n_wells) in Mt/year
- `optimal_storage::Float64`: Total CO2 stored at end of simulation (Mt)
- `optimal_leakage::Float64`: Total CO2 leaked at end of simulation (Mt)
- `time_periods::Vector{Float64}`: Start time of each period (years)
- `objective_value::Float64`: Final objective function value
- `convergence_info::Dict{Symbol, Any}`: Algorithm-specific convergence information
- `evaluation_history::Vector{Tuple{Vector{Float64}, Float64}}`: History of (params, objective)
"""
struct InjectionOptimizationResult
    optimal_rates::Matrix{Float64}
    optimal_storage::Float64
    optimal_leakage::Float64
    time_periods::Vector{Float64}
    objective_value::Float64
    convergence_info::Dict{Symbol, Any}
    evaluation_history::Vector{Tuple{Vector{Float64}, Float64}}
end


# =============================================================================
# MULTI-WELL INJECTION EVENT GENERATION
# =============================================================================

"""
    create_multi_well_injection_events(
        rates::Matrix{Float64},
        config::OptimizationConfig,
        layers::Vector{Layer}
    ) -> Vector{Vector{InjectionEvent}}

Convert optimization decision variables (rates matrix) to injection events for each layer.

# Arguments
- `rates`: Matrix of injection rates (n_periods x n_wells) in Mt/year
- `config`: OptimizationConfig with well locations and time settings
- `layers`: Vector of Layer structs for the reservoir

# Returns
Vector of InjectionEvent vectors, one per layer. Only layers with wells have non-zero injection.
"""
function create_multi_well_injection_events(
    rates::Matrix{Float64},
    config::OptimizationConfig,
    layers::Vector{Layer}
)::Vector{Vector{InjectionEvent}}

    n_periods = config.n_time_periods
    n_layers = length(layers)

    # Time period boundaries
    period_duration = (config.end_time - config.start_time) / n_periods
    time_starts = [config.start_time + (i - 1) * period_duration for i in 1:n_periods]

    # Grid size from bottom layer
    grid_size = size(layers[1].trap_structure.topography)

    # Initialize injection events for each layer with zero injection
    zero_injection = zeros(grid_size)
    injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)
    for layer_idx in 1:n_layers
        injection_events[layer_idx] = [InjectionEvent(0.0, zero_injection)]
    end

    # Group wells by layer
    wells_by_layer = Dict{Int, Vector{Int}}()
    for (well_idx, layer_idx) in enumerate(config.well_layers)
        if !haskey(wells_by_layer, layer_idx)
            wells_by_layer[layer_idx] = Int[]
        end
        push!(wells_by_layer[layer_idx], well_idx)
    end

    # Create injection events for each layer with wells
    for (layer_idx, well_indices) in wells_by_layer
        layer_events = Vector{InjectionEvent}(undef, n_periods)

        for period_idx in 1:n_periods
            # Create spatial injection matrix
            injection_rate_matrix = zeros(grid_size)

            for well_idx in well_indices
                # Convert Mt/year to m³/year
                rate_mt_per_year = rates[period_idx, well_idx]
                rate_m3_per_year = rate_mt_per_year * 1e9 / config.co2_density

                # Place injection at well location
                well_loc = config.well_locations[well_idx]
                injection_rate_matrix[well_loc] = rate_m3_per_year
            end

            layer_events[period_idx] = InjectionEvent(time_starts[period_idx], injection_rate_matrix)
        end

        injection_events[layer_idx] = layer_events
    end

    return injection_events
end


# =============================================================================
# OBJECTIVE FUNCTION AND EVALUATION
# =============================================================================

"""
    evaluate_injection_schedule(
        rates::Matrix{Float64},
        config::OptimizationConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
        return_snapshots::Bool = false
    )

Evaluate an injection schedule by running the full simulation.

# Returns
- `storage`: Total CO2 stored at end of simulation (Mt)
- `leakage`: Total CO2 leaked at end of simulation (Mt)
- `snapshots`: (optional) Vector of ReservoirSnapshot for analysis
"""
function evaluate_injection_schedule(
    rates::Matrix{Float64},
    config::OptimizationConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
    return_snapshots::Bool = false
)
    # Create injection events from rates
    injection_events = create_multi_well_injection_events(rates, config, layers)

    # Run simulation
    seqs, leakage_states = fill_layers(
        layers, domain, reservoir_properties, injection_events;
        verbose=false
    )

    # Generate snapshots (just start and end for efficiency)
    snapshots = generate_reservoir_snapshots(
        layers, seqs, leakage_states, domain, reservoir_properties, injection_events;
        num_snapshots=2,
        start_time=config.start_time,
        end_time=config.end_time,
        verbose=false
    )

    # Extract final state
    final_snapshot = snapshots[end]

    # Convert m³ to Mt
    storage_mt = final_snapshot.total_stored_m3 * config.co2_density / 1e9
    leakage_mt = final_snapshot.total_leaked_m3 * config.co2_density / 1e9

    if return_snapshots
        return (storage_mt, leakage_mt, snapshots)
    else
        return (storage_mt, leakage_mt)
    end
end


"""
    create_objective_function(
        config::OptimizationConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
        mass_penalty_weight::Float64 = 1000.0
    ) -> Function

Create the objective function for optimization.

The objective function:
1. Converts flat decision vector to rates matrix
2. Evaluates the injection schedule
3. Applies penalty for mass constraint violation
4. Returns negative storage (for minimization)
"""
function create_objective_function(
    config::OptimizationConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
    mass_penalty_weight::Float64 = 1000.0
)
    n_periods = config.n_time_periods
    n_wells = length(config.well_locations)
    period_duration = (config.end_time - config.start_time) / n_periods

    function objective(x::Vector{Float64})
        # Reshape to rates matrix
        rates = reshape(x, n_periods, n_wells)

        # Compute total injected mass
        total_mass = sum(rates) * period_duration

        # Mass constraint penalty (squared error)
        mass_error = abs(total_mass - config.total_mass_mt)
        mass_penalty = mass_penalty_weight * mass_error^2

        # Run simulation
        try
            storage, _ = evaluate_injection_schedule(
                rates, config, layers, domain, reservoir_properties
            )

            # Objective: minimize negative storage + penalty
            return -storage + mass_penalty

        catch e
            # If simulation fails, return large penalty
            @warn "Simulation failed: $e"
            return 1e10
        end
    end

    return objective
end


# =============================================================================
# OPTIMIZATION BACKENDS
# =============================================================================

"""
Run optimization using BlackBoxOptim.jl (Differential Evolution).
"""
function _run_blackbox_optimization(
    objective::Function,
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64},
    method::Symbol,
    max_iterations::Int,
    verbose::Bool
)
    # Create search range
    search_range = [(lower_bounds[i], upper_bounds[i]) for i in eachindex(lower_bounds)]

    # Run optimization
    result = bboptimize(
        objective;
        SearchRange = search_range,
        Method = method,
        MaxSteps = max_iterations,
        TraceMode = verbose ? :verbose : :silent
    )

    return Dict(
        :optimal_x => best_candidate(result),
        :optimal_obj => best_fitness(result),
        :convergence_info => Dict(
            :algorithm => method,
            :n_evaluations => max_iterations
        )
    )
end


"""
Simple random search for baseline comparison.
"""
function _run_random_search(
    objective::Function,
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64},
    max_iterations::Int,
    verbose::Bool
)
    n_vars = length(lower_bounds)
    best_x = lower_bounds + rand(n_vars) .* (upper_bounds - lower_bounds)
    best_obj = objective(best_x)

    for iter in 1:max_iterations
        x = lower_bounds + rand(n_vars) .* (upper_bounds - lower_bounds)
        obj = objective(x)

        if obj < best_obj
            best_x = x
            best_obj = obj
            verbose && println("Iteration $iter: new best = $(round(best_obj, digits=4))")
        end
    end

    return Dict(
        :optimal_x => best_x,
        :optimal_obj => best_obj,
        :convergence_info => Dict(
            :algorithm => :random_search,
            :n_evaluations => max_iterations
        )
    )
end


"""
Latin Hypercube sampling for initial points.
"""
function _latin_hypercube_sample(n_samples::Int, lower_bounds::Vector{Float64}, upper_bounds::Vector{Float64})
    n_vars = length(lower_bounds)
    samples = zeros(n_samples, n_vars)

    for j in 1:n_vars
        # Create evenly spaced intervals
        intervals = [(i-1)/n_samples + rand()/n_samples for i in 1:n_samples]
        shuffle!(intervals)

        # Scale to bounds
        for i in 1:n_samples
            samples[i, j] = lower_bounds[j] + intervals[i] * (upper_bounds[j] - lower_bounds[j])
        end
    end

    return samples
end


"""
Run Bayesian-like optimization using surrogate model approach.
Uses Latin Hypercube sampling + local refinement around best points.
"""
function _run_bayesian_optimization(
    objective::Function,
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64},
    n_vars::Int,
    config::BayesianOptimizationConfig,
    verbose::Bool
)
    # Initial Latin Hypercube sampling
    verbose && println("Generating $(config.n_initial_samples) initial samples...")
    initial_samples = _latin_hypercube_sample(config.n_initial_samples, lower_bounds, upper_bounds)

    # Evaluate initial samples
    sample_values = Float64[]
    sample_points = Vector{Vector{Float64}}()

    for i in 1:config.n_initial_samples
        x = initial_samples[i, :]
        obj_val = objective(x)
        push!(sample_points, x)
        push!(sample_values, obj_val)
        verbose && println("Initial sample $i/$(config.n_initial_samples): obj = $(round(obj_val, digits=4))")
    end

    best_idx = argmin(sample_values)
    best_x = sample_points[best_idx]
    best_obj = sample_values[best_idx]

    verbose && println("\nBest initial: $(round(best_obj, digits=4))")

    # Iterative refinement around best points
    n_remaining = config.max_evaluations - config.n_initial_samples
    no_improvement_count = 0
    exploration_radius = 0.2  # Start with 20% of range

    for iter in 1:n_remaining
        # Adaptive exploration: shrink radius as we converge
        current_radius = exploration_radius * (1.0 - 0.5 * iter / n_remaining)

        # Sample near the best point (exploitation) or random (exploration)
        if rand() < 0.7  # 70% exploitation, 30% exploration
            # Sample near best point
            perturbation = randn(n_vars) .* current_radius .* (upper_bounds - lower_bounds)
            x_new = clamp.(best_x + perturbation, lower_bounds, upper_bounds)
        else
            # Random exploration
            x_new = lower_bounds + rand(n_vars) .* (upper_bounds - lower_bounds)
        end

        obj_new = objective(x_new)
        push!(sample_points, x_new)
        push!(sample_values, obj_new)

        verbose && println("Iteration $iter: obj = $(round(obj_new, digits=4)), best = $(round(best_obj, digits=4))")

        # Check for improvement
        if obj_new < best_obj - config.convergence_threshold
            best_x = x_new
            best_obj = obj_new
            no_improvement_count = 0
        else
            no_improvement_count += 1
        end

        # Early stopping
        if no_improvement_count >= config.convergence_patience
            verbose && println("Converged after $iter iterations (no improvement)")
            break
        end
    end

    return Dict(
        :optimal_x => collect(best_x),
        :optimal_obj => best_obj,
        :convergence_info => Dict(
            :algorithm => :bayesian,
            :n_evaluations => length(sample_values),
            :converged => no_improvement_count >= config.convergence_patience
        )
    )
end


# =============================================================================
# MAIN OPTIMIZATION INTERFACE
# =============================================================================

"""
    optimize_injection_schedule(
        config::OptimizationConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
        algorithm::Symbol = :bayesian,
        bo_config::Union{BayesianOptimizationConfig, Nothing} = nothing,
        max_iterations::Int = 100,
        verbose::Bool = true,
        seed::Union{Int, Nothing} = nothing
    ) -> InjectionOptimizationResult

Optimize injection schedule to maximize CO2 storage.

# Arguments
- `config`: OptimizationConfig defining the problem
- `layers`: Reservoir layer structure
- `domain`: Spatial domain
- `reservoir_properties`: Reservoir properties
- `algorithm`: Optimization algorithm (:bayesian, :differential_evolution, :random_search)
- `bo_config`: Bayesian optimization configuration (for :bayesian algorithm)
- `max_iterations`: Maximum optimization iterations
- `verbose`: Print progress
- `seed`: Random seed for reproducibility

# Returns
InjectionOptimizationResult with optimal schedule and metadata
"""
function optimize_injection_schedule(
    config::OptimizationConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
    algorithm::Symbol = :bayesian,
    bo_config::Union{BayesianOptimizationConfig, Nothing} = nothing,
    max_iterations::Int = 100,
    verbose::Bool = true,
    seed::Union{Int, Nothing} = nothing
)
    # Set random seed
    if !isnothing(seed)
        Random.seed!(seed)
    end

    n_periods = config.n_time_periods
    n_wells = length(config.well_locations)
    n_vars = n_periods * n_wells

    # Variable bounds
    lower_bounds = fill(config.min_rate_mt_per_year, n_vars)
    upper_bounds = fill(config.max_rate_mt_per_year, n_vars)

    # Create objective function
    objective = create_objective_function(config, layers, domain, reservoir_properties)

    # Track evaluation history
    evaluation_history = Vector{Tuple{Vector{Float64}, Float64}}()

    # Wrapper to track history
    function tracked_objective(x)
        obj_val = objective(x)
        push!(evaluation_history, (copy(x), obj_val))
        return obj_val
    end

    verbose && println("\nStarting optimization with $algorithm algorithm")
    verbose && println("Decision variables: $n_vars ($n_periods periods × $n_wells wells)")
    verbose && println("Target mass: $(config.total_mass_mt) Mt")
    verbose && println("")

    # Run optimization based on algorithm choice
    if algorithm == :bayesian
        result = _run_bayesian_optimization(
            tracked_objective, lower_bounds, upper_bounds, n_vars,
            isnothing(bo_config) ? BayesianOptimizationConfig(max_evaluations=max_iterations) : bo_config,
            verbose
        )
    elseif algorithm == :differential_evolution || algorithm == :adaptive_de_rand_1_bin
        result = _run_blackbox_optimization(
            tracked_objective, lower_bounds, upper_bounds,
            :adaptive_de_rand_1_bin, max_iterations, verbose
        )
    elseif algorithm == :de_rand_1_bin
        result = _run_blackbox_optimization(
            tracked_objective, lower_bounds, upper_bounds,
            :de_rand_1_bin, max_iterations, verbose
        )
    elseif algorithm == :random_search
        result = _run_random_search(
            tracked_objective, lower_bounds, upper_bounds,
            max_iterations, verbose
        )
    else
        error("Unknown algorithm: $algorithm. Options: :bayesian, :differential_evolution, :de_rand_1_bin, :random_search")
    end

    # Extract optimal solution
    optimal_x = result[:optimal_x]
    optimal_rates = reshape(optimal_x, n_periods, n_wells)

    # Evaluate final solution for accurate metrics
    storage, leakage = evaluate_injection_schedule(
        optimal_rates, config, layers, domain, reservoir_properties
    )

    # Time periods
    period_duration = (config.end_time - config.start_time) / n_periods
    time_periods = [config.start_time + (i - 1) * period_duration for i in 1:n_periods]

    verbose && println("\nOptimization complete!")
    verbose && println("  Storage: $(round(storage, digits=3)) Mt")
    verbose && println("  Leakage: $(round(leakage, digits=3)) Mt")
    verbose && println("  Efficiency: $(round(100 * storage / config.total_mass_mt, digits=1))%")

    return InjectionOptimizationResult(
        optimal_rates,
        storage,
        leakage,
        time_periods,
        result[:optimal_obj],
        result[:convergence_info],
        evaluation_history
    )
end


# =============================================================================
# CONVENIENCE FUNCTION FOR SLEIPNER
# =============================================================================

"""
    optimize_sleipner_injection(;
        n_time_periods::Int = 10,
        total_mass_mt::Float64 = 12.18,
        well2_location::Union{CartesianIndex, Nothing} = nothing,
        well2_layer::Int = 1,
        max_rate_mt_per_year::Float64 = 2.0,
        algorithm::Symbol = :bayesian,
        max_evaluations::Int = 100,
        verbose::Bool = true,
        seed::Union{Int, Nothing} = nothing
    ) -> InjectionOptimizationResult

Convenience function to run injection optimization for Sleipner with sensible defaults.

Sets up the full reservoir, defines two wells, and runs optimization.

# Arguments
- `n_time_periods`: Number of time periods (default: 10)
- `total_mass_mt`: Total mass to inject in Mt (default: 12.18 = Sleipner historical)
- `well2_location`: Grid location for second well (default: offset from main feeder)
- `well2_layer`: Layer for second well (default: 1)
- `max_rate_mt_per_year`: Maximum injection rate per well (Mt/year)
- `algorithm`: Optimization algorithm (:bayesian, :differential_evolution)
- `max_evaluations`: Maximum objective evaluations
- `verbose`: Print progress
- `seed`: Random seed for reproducibility

# Example
```julia
result = optimize_sleipner_injection(
    n_time_periods = 10,
    well2_layer = 3,  # Second well in layer 3
    algorithm = :bayesian,
    max_evaluations = 50
)
println("Storage: \$(result.optimal_storage) Mt")
```
"""
function optimize_sleipner_injection(;
    n_time_periods::Int = 10,
    total_mass_mt::Float64 = 12.2,
    well2_location::Union{CartesianIndex, Nothing} = nothing,
    well2_layer::Int = 1,
    max_rate_mt_per_year::Float64 = 2.0,
    algorithm::Symbol = :bayesian,
    max_evaluations::Int = 100,
    verbose::Bool = true,
    seed::Union{Int, Nothing} = nothing
)
    # Load reservoir
    verbose && println("Loading Sleipner topography...")
    topography = load_sleipner_topography()
    domain = create_domain_from_topography(topography, 1.0)
    layers = analyze_base_surfaces(topography; boundary_condition=:closed)
    reservoir_properties = generate_reservoir_properties_for_sleipner_layers()

    # Set up wells
    well1, (utm_x1, utm_y1, _) = load_feeder_location(topography)

    # Default second well: offset 500m east from main feeder
    if isnothing(well2_location)
        well2 = utm_to_grid_index(utm_x1 + 500.0, utm_y1, topography)
    else
        well2 = well2_location
    end

    verbose && println("Well 1: $well1 (layer 1, main feeder)")
    verbose && println("Well 2: $well2 (layer $well2_layer)")

    # Create config
    config = OptimizationConfig(
        n_time_periods,
        [well1, well2],
        [1, well2_layer],
        total_mass_mt,
        max_rate_mt_per_year,
        0.0,  # min rate
        0.0,  # start time
        15.0, # end time
        reservoir_properties[1].co2_density # CO2 density
    )

    # Run optimization
    result = optimize_injection_schedule(
        config, layers, domain, reservoir_properties;
        algorithm = algorithm,
        bo_config = algorithm == :bayesian ? BayesianOptimizationConfig(max_evaluations=max_evaluations) : nothing,
        max_iterations = max_evaluations,
        verbose = verbose,
        seed = seed
    )

    return result
end


# =============================================================================
# VISUALIZATION
# =============================================================================

"""
    plot_injection_schedule(
        result::InjectionOptimizationResult,
        config::OptimizationConfig;
        output_file::Union{String, Nothing} = nothing,
        figsize::Tuple{Int, Int} = (1000, 600)
    )

Plot the optimal injection schedule as a grouped bar chart.
"""
function plot_injection_schedule(
    result::InjectionOptimizationResult,
    config::OptimizationConfig;
    output_file::Union{String, Nothing} = nothing,
    figsize::Tuple{Int, Int} = (1000, 600)
)
    n_periods = config.n_time_periods
    n_wells = length(config.well_locations)

    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1],
        xlabel = "Time Period Start (years)",
        ylabel = "Injection Rate (Mt/year)",
        title = "Optimal Injection Schedule"
    )

    # Colors for wells
    colors = [:blue, :red, :green, :orange, :purple]

    # Bar width and positions
    bar_width = 0.8 / n_wells
    period_duration = (config.end_time - config.start_time) / n_periods

    for well_idx in 1:n_wells
        # Center bars for each period
        offset = (well_idx - (n_wells + 1) / 2) * bar_width
        positions = result.time_periods .+ offset

        rates = result.optimal_rates[:, well_idx]
        layer = config.well_layers[well_idx]
        barplot!(ax, positions, rates, width=bar_width * period_duration * 0.9,
                 color=colors[mod1(well_idx, length(colors))],
                 label="Well $well_idx (L$layer)")
    end

    axislegend(ax, position=:rt)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Saved: $output_file")
    end

    return fig
end


"""
    plot_optimization_convergence(
        result::InjectionOptimizationResult;
        output_file::Union{String, Nothing} = nothing,
        figsize::Tuple{Int, Int} = (800, 500)
    )

Plot the convergence history of the optimization.
"""
function plot_optimization_convergence(
    result::InjectionOptimizationResult;
    output_file::Union{String, Nothing} = nothing,
    figsize::Tuple{Int, Int} = (800, 500)
)
    history = result.evaluation_history
    n_evals = length(history)

    objectives = [h[2] for h in history]
    best_so_far = accumulate(min, objectives)

    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1],
        xlabel = "Evaluation",
        ylabel = "Objective Value (negative storage + penalty)",
        title = "Optimization Convergence"
    )

    scatter!(ax, 1:n_evals, objectives, color=(:gray, 0.5), markersize=5, label="Evaluations")
    lines!(ax, 1:n_evals, best_so_far, color=:blue, linewidth=2, label="Best so far")

    axislegend(ax, position=:rt)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Saved: $output_file")
    end

    return fig
end


# =============================================================================
# WELL LOCATION OPTIMIZATION
# =============================================================================

"""
    LocationOptimizationConfig

Configuration for optimizing well placement locations.

# Fields
- `n_wells::Int`: Number of wells to place
- `n_layers::Int`: Number of layers in the reservoir (typically 9 for Sleipner)
- `grid_size::Tuple{Int, Int}`: Grid dimensions (nx, ny)
- `allowed_layers::Vector{Int}`: Which layers wells can be placed in (default: all)
- `total_mass_mt::Float64`: Total CO2 mass to inject (Mt)
- `injection_rate_mt_per_year::Float64`: Injection rate per well (Mt/year) - uniform for location optimization
- `start_time::Float64`: Simulation start time (years)
- `end_time::Float64`: Simulation end time (years)
- `co2_density::Float64`: CO2 density for volume conversion (kg/m³)
- `min_well_spacing::Int`: Minimum grid cells between wells (default: 5)
- `require_bottom_layer_well::Bool`: If true, force at least one well to be in layer 1 (default: false)
"""
struct LocationOptimizationConfig
    n_wells::Int
    n_layers::Int
    grid_size::Tuple{Int, Int}
    allowed_layers::Vector{Int}
    total_mass_mt::Float64
    injection_rate_mt_per_year::Float64
    start_time::Float64
    end_time::Float64
    co2_density::Float64
    min_well_spacing::Int
    require_bottom_layer_well::Bool

    function LocationOptimizationConfig(;
        n_wells::Int,
        n_layers::Int = 9,
        grid_size::Tuple{Int, Int},
        allowed_layers::Union{Vector{Int}, Nothing} = nothing,
        total_mass_mt::Float64 = 12.18,
        injection_rate_mt_per_year::Union{Float64, Nothing} = nothing,
        start_time::Float64 = 0.0,
        end_time::Float64 = 15.0,
        co2_density::Float64 = 425.0,
        min_well_spacing::Int = 5,
        require_bottom_layer_well::Bool = false
    )
        @assert n_wells >= 1 "Need at least 1 well"
        @assert n_layers >= 1 "Need at least 1 layer"
        @assert grid_size[1] >= 1 && grid_size[2] >= 1 "Grid size must be positive"

        # Default: all layers except top (caprock)
        if isnothing(allowed_layers)
            allowed_layers = collect(1:(n_layers-1))
        end

        # If requiring bottom layer well, ensure layer 1 is in allowed_layers
        if require_bottom_layer_well && !(1 in allowed_layers)
            @warn "require_bottom_layer_well=true but layer 1 not in allowed_layers. Adding layer 1."
            allowed_layers = sort(unique([1; allowed_layers]))
        end

        # Default injection rate: distribute total mass evenly
        if isnothing(injection_rate_mt_per_year)
            injection_rate_mt_per_year = total_mass_mt / (n_wells * (end_time - start_time))
        end

        new(n_wells, n_layers, grid_size, allowed_layers, total_mass_mt,
            injection_rate_mt_per_year, start_time, end_time, co2_density, min_well_spacing,
            require_bottom_layer_well)
    end
end


"""
    LocationOptimizationResult

Result from well location optimization.

# Fields
- `optimal_locations::Vector{CartesianIndex{2}}`: Optimal grid locations for each well
- `optimal_layers::Vector{Int}`: Optimal layer for each well
- `optimal_storage::Float64`: Total CO2 stored at end of simulation (Mt)
- `optimal_leakage::Float64`: Total CO2 leaked at end of simulation (Mt)
- `storage_efficiency::Float64`: Fraction of injected CO2 that is stored
- `objective_value::Float64`: Final objective function value
- `convergence_info::Dict{Symbol, Any}`: Algorithm-specific convergence information
- `evaluation_history::Vector{Tuple{Vector{Float64}, Float64}}`: History of (params, objective)
"""
struct LocationOptimizationResult
    optimal_locations::Vector{CartesianIndex{2}}
    optimal_layers::Vector{Int}
    optimal_storage::Float64
    optimal_leakage::Float64
    storage_efficiency::Float64
    objective_value::Float64
    convergence_info::Dict{Symbol, Any}
    evaluation_history::Vector{Tuple{Vector{Float64}, Float64}}
end


"""
    decode_location_variables(x::Vector{Float64}, config::LocationOptimizationConfig)

Convert normalized decision variables [0,1] to well locations and layers.

Each well uses 3 variables: (x_norm, y_norm, layer_norm)

If `config.require_bottom_layer_well` is true, the first well is forced to layer 1.
"""
function decode_location_variables(x::Vector{Float64}, config::LocationOptimizationConfig)
    n_wells = config.n_wells
    nx, ny = config.grid_size

    locations = Vector{CartesianIndex{2}}(undef, n_wells)
    layers = Vector{Int}(undef, n_wells)

    for w in 1:n_wells
        idx = (w - 1) * 3
        x_norm = x[idx + 1]
        y_norm = x[idx + 2]
        layer_norm = x[idx + 3]

        # Map to grid indices (1-based, with margin from edges)
        margin = 2
        i = clamp(round(Int, margin + x_norm * (nx - 2*margin)), 1, nx)
        j = clamp(round(Int, margin + y_norm * (ny - 2*margin)), 1, ny)
        locations[w] = CartesianIndex(i, j)

        # Map to allowed layers
        # If require_bottom_layer_well is true, force the first well to layer 1
        if config.require_bottom_layer_well && w == 1
            layers[w] = 1
        else
            n_allowed = length(config.allowed_layers)
            layer_idx = clamp(round(Int, 1 + layer_norm * (n_allowed - 1)), 1, n_allowed)
            layers[w] = config.allowed_layers[layer_idx]
        end
    end

    return locations, layers
end


"""
    compute_well_spacing_penalty(locations::Vector{CartesianIndex{2}}, min_spacing::Int)

Compute penalty for wells that are too close together.
"""
function compute_well_spacing_penalty(locations::Vector{CartesianIndex{2}}, min_spacing::Int)
    penalty = 0.0
    n_wells = length(locations)

    for i in 1:n_wells
        for j in (i+1):n_wells
            dist = sqrt((locations[i][1] - locations[j][1])^2 +
                       (locations[i][2] - locations[j][2])^2)
            if dist < min_spacing
                penalty += (min_spacing - dist)^2
            end
        end
    end

    return penalty
end


"""
    evaluate_well_locations(
        locations::Vector{CartesianIndex{2}},
        layers::Vector{Int},
        config::LocationOptimizationConfig,
        reservoir_layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
        return_snapshots::Bool = false
    )

Evaluate a set of well locations by running the simulation with uniform injection rates.
"""
function evaluate_well_locations(
    locations::Vector{CartesianIndex{2}},
    layers::Vector{Int},
    config::LocationOptimizationConfig,
    reservoir_layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
    return_snapshots::Bool = false
)
    n_layers = length(reservoir_layers)
    grid_size = size(reservoir_layers[1].trap_structure.topography)

    # Create uniform injection rate over simulation period
    rate_per_well = config.injection_rate_mt_per_year
    rate_m3_per_year = rate_per_well * 1e9 / config.co2_density

    # Group wells by layer
    wells_by_layer = Dict{Int, Vector{Int}}()
    for (well_idx, layer_idx) in enumerate(layers)
        if !haskey(wells_by_layer, layer_idx)
            wells_by_layer[layer_idx] = Int[]
        end
        push!(wells_by_layer[layer_idx], well_idx)
    end

    # Create injection events for each layer
    zero_injection = zeros(grid_size)
    injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)
    for layer_idx in 1:n_layers
        injection_events[layer_idx] = [InjectionEvent(0.0, zero_injection)]
    end

    # Add injection for layers with wells
    for (layer_idx, well_indices) in wells_by_layer
        injection_rate_matrix = zeros(grid_size)
        for well_idx in well_indices
            loc = locations[well_idx]
            injection_rate_matrix[loc] = rate_m3_per_year
        end
        injection_events[layer_idx] = [InjectionEvent(config.start_time, injection_rate_matrix)]
    end

    # Run simulation
    seqs, leakage_states = fill_layers(
        reservoir_layers, domain, reservoir_properties, injection_events;
        verbose=false
    )

    # Generate snapshots
    snapshots = generate_reservoir_snapshots(
        reservoir_layers, seqs, leakage_states, domain, reservoir_properties, injection_events;
        num_snapshots=2,
        start_time=config.start_time,
        end_time=config.end_time,
        verbose=false
    )

    # Extract final state
    final_snapshot = snapshots[end]

    # Convert m³ to Mt
    storage_mt = final_snapshot.total_stored_m3 * config.co2_density / 1e9
    leakage_mt = final_snapshot.total_leaked_m3 * config.co2_density / 1e9

    if return_snapshots
        return (storage_mt, leakage_mt, snapshots)
    else
        return (storage_mt, leakage_mt)
    end
end


"""
    create_location_objective_function(
        config::LocationOptimizationConfig,
        reservoir_layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
        spacing_penalty_weight::Float64 = 100.0
    ) -> Function

Create the objective function for well location optimization.
"""
function create_location_objective_function(
    config::LocationOptimizationConfig,
    reservoir_layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
    spacing_penalty_weight::Float64 = 100.0
)
    function objective(x::Vector{Float64})
        # Decode decision variables to locations and layers
        locations, layers = decode_location_variables(x, config)

        # Compute spacing penalty
        spacing_penalty = spacing_penalty_weight * compute_well_spacing_penalty(locations, config.min_well_spacing)

        # Run simulation
        try
            storage, _ = evaluate_well_locations(
                locations, layers, config, reservoir_layers, domain, reservoir_properties
            )

            # Objective: minimize negative storage + penalties
            return -storage + spacing_penalty

        catch e
            @warn "Simulation failed: $e"
            return 1e10
        end
    end

    return objective
end


"""
    optimize_well_locations(
        config::LocationOptimizationConfig,
        layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
        algorithm::Symbol = :differential_evolution,
        max_iterations::Int = 100,
        verbose::Bool = true,
        seed::Union{Int, Nothing} = nothing
    ) -> LocationOptimizationResult

Optimize well placement locations to maximize CO2 storage.

# Arguments
- `config`: LocationOptimizationConfig defining the problem
- `layers`: Reservoir layer structure
- `domain`: Spatial domain
- `reservoir_properties`: Reservoir properties
- `algorithm`: Optimization algorithm (:bayesian, :differential_evolution, :random_search)
- `max_iterations`: Maximum optimization iterations
- `verbose`: Print progress
- `seed`: Random seed for reproducibility

# Returns
LocationOptimizationResult with optimal locations and metadata

# Example
```julia
config = LocationOptimizationConfig(
    n_wells = 2,
    grid_size = (64, 118),
    allowed_layers = [1, 2, 3],
    total_mass_mt = 12.18
)
result = optimize_well_locations(config, layers, domain, reservoir_properties)
```
"""
function optimize_well_locations(
    config::LocationOptimizationConfig,
    layers::Vector{Layer},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}};
    algorithm::Symbol = :differential_evolution,
    bo_config::Union{BayesianOptimizationConfig, Nothing} = nothing,
    max_iterations::Int = 100,
    verbose::Bool = true,
    seed::Union{Int, Nothing} = nothing
)
    # Set random seed
    if !isnothing(seed)
        Random.seed!(seed)
    end

    n_vars = config.n_wells * 3  # 3 variables per well: x, y, layer

    # All variables normalized to [0, 1]
    lower_bounds = zeros(n_vars)
    upper_bounds = ones(n_vars)

    # Create objective function
    objective = create_location_objective_function(config, layers, domain, reservoir_properties)

    # Track evaluation history
    evaluation_history = Vector{Tuple{Vector{Float64}, Float64}}()

    function tracked_objective(x)
        obj_val = objective(x)
        push!(evaluation_history, (copy(x), obj_val))
        return obj_val
    end

    verbose && println("\nStarting well location optimization")
    verbose && println("Number of wells: $(config.n_wells)")
    verbose && println("Allowed layers: $(config.allowed_layers)")
    if config.require_bottom_layer_well
        verbose && println("Constraint: At least one well must be in layer 1 (bottom)")
    end
    verbose && println("Grid size: $(config.grid_size)")
    verbose && println("Decision variables: $n_vars (3 per well: x, y, layer)")
    verbose && println("")

    # Run optimization
    if algorithm == :bayesian
        result = _run_bayesian_optimization(
            tracked_objective, lower_bounds, upper_bounds, n_vars,
            isnothing(bo_config) ? BayesianOptimizationConfig(max_evaluations=max_iterations) : bo_config,
            verbose
        )
    elseif algorithm == :differential_evolution || algorithm == :adaptive_de_rand_1_bin
        result = _run_blackbox_optimization(
            tracked_objective, lower_bounds, upper_bounds,
            :adaptive_de_rand_1_bin, max_iterations, verbose
        )
    elseif algorithm == :random_search
        result = _run_random_search(
            tracked_objective, lower_bounds, upper_bounds,
            max_iterations, verbose
        )
    else
        error("Unknown algorithm: $algorithm")
    end

    # Decode optimal solution
    optimal_x = result[:optimal_x]
    optimal_locations, optimal_layers = decode_location_variables(optimal_x, config)

    # Evaluate final solution
    storage, leakage = evaluate_well_locations(
        optimal_locations, optimal_layers, config, layers, domain, reservoir_properties
    )

    total_injected = config.n_wells * config.injection_rate_mt_per_year * (config.end_time - config.start_time)
    efficiency = storage / total_injected

    verbose && println("\nLocation optimization complete!")
    verbose && println("  Optimal locations:")
    for (i, (loc, layer)) in enumerate(zip(optimal_locations, optimal_layers))
        verbose && println("    Well $i: $loc in Layer $layer")
    end
    verbose && println("  Storage: $(round(storage, digits=3)) Mt")
    verbose && println("  Leakage: $(round(leakage, digits=3)) Mt")
    verbose && println("  Efficiency: $(round(100 * efficiency, digits=1))%")

    return LocationOptimizationResult(
        optimal_locations,
        optimal_layers,
        storage,
        leakage,
        efficiency,
        result[:optimal_obj],
        result[:convergence_info],
        evaluation_history
    )
end


"""
    optimize_sleipner_well_locations(;
        n_wells::Int = 2,
        allowed_layers::Union{Vector{Int}, Nothing} = nothing,
        require_bottom_layer_well::Bool = false,
        total_mass_mt::Float64 = 12.18,
        algorithm::Symbol = :differential_evolution,
        max_evaluations::Int = 100,
        verbose::Bool = true,
        seed::Union{Int, Nothing} = nothing
    ) -> LocationOptimizationResult

Convenience function to optimize well locations for Sleipner.

# Arguments
- `n_wells`: Number of wells to place (default: 2)
- `allowed_layers`: Which layers to consider (default: 1-8, excluding caprock)
- `require_bottom_layer_well`: If true, force at least one well to be in layer 1 (default: false)
- `total_mass_mt`: Total mass to inject in Mt (default: 12.18)
- `algorithm`: Optimization algorithm (:differential_evolution recommended for this problem)
- `max_evaluations`: Maximum objective evaluations
- `verbose`: Print progress
- `seed`: Random seed for reproducibility

# Example
```julia
result = optimize_sleipner_well_locations(
    n_wells = 3,
    allowed_layers = [1, 2, 3],  # Only lower layers
    require_bottom_layer_well = true,  # Ensure one well in L1
    max_evaluations = 50
)
```
"""
function optimize_sleipner_well_locations(;
    n_wells::Int = 2,
    allowed_layers::Union{Vector{Int}, Nothing} = nothing,
    require_bottom_layer_well::Bool = false,
    total_mass_mt::Float64 = 12.2,
    algorithm::Symbol = :differential_evolution,
    max_evaluations::Int = 100,
    verbose::Bool = true,
    seed::Union{Int, Nothing} = nothing
)
    # Load reservoir
    verbose && println("Loading Sleipner topography...")
    topography = load_sleipner_topography()
    domain = create_domain_from_topography(topography, 1.0)
    layers = analyze_base_surfaces(topography; boundary_condition=:closed)
    reservoir_properties = generate_reservoir_properties_for_sleipner_layers()

    # Create config
    config = LocationOptimizationConfig(
        n_wells = n_wells,
        n_layers = length(layers),
        grid_size = (topography.nx, topography.ny),
        allowed_layers = allowed_layers,
        total_mass_mt = total_mass_mt,
        start_time = 0.0,
        end_time = 15.0,
        co2_density = 460.0,
        min_well_spacing = 5,
        require_bottom_layer_well = require_bottom_layer_well
    )

    # Run optimization
    result = optimize_well_locations(
        config, layers, domain, reservoir_properties;
        algorithm = algorithm,
        max_iterations = max_evaluations,
        verbose = verbose,
        seed = seed
    )

    return result
end


"""
    plot_well_locations(
        result::LocationOptimizationResult,
        layers::Vector{Layer};
        output_file::Union{String, Nothing} = nothing,
        figsize::Tuple{Int, Int} = (1000, 800),
        show_topography::Bool = true
    )

Plot the optimal well locations on the reservoir grid.
"""
function plot_well_locations(
    result::LocationOptimizationResult,
    layers::Vector{Layer};
    output_file::Union{String, Nothing} = nothing,
    figsize::Tuple{Int, Int} = (1000, 800),
    show_topography::Bool = true
)
    fig = Figure(size=figsize)

    # Get topography from bottom layer for background
    topo = layers[1].trap_structure.topography

    ax = Axis(fig[1, 1],
        xlabel = "X (grid index)",
        ylabel = "Y (grid index)",
        title = "Optimal Well Locations\nStorage: $(round(result.optimal_storage, digits=2)) Mt, Efficiency: $(round(100*result.storage_efficiency, digits=1))%",
        aspect = DataAspect()
    )

    if show_topography
        hm = heatmap!(ax, topo, colormap=:viridis)
        Colorbar(fig[1, 2], hm, label="Depth (m)")
    end

    # Plot well locations
    colors = [:red, :orange, :yellow, :cyan, :magenta]
    for (i, (loc, layer)) in enumerate(zip(result.optimal_locations, result.optimal_layers))
        scatter!(ax, [loc[1]], [loc[2]],
                color=colors[mod1(i, length(colors))],
                markersize=20,
                marker=:star5,
                strokewidth=2,
                strokecolor=:white,
                label="Well $i (L$layer)")
    end

    axislegend(ax, position=:rt)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Saved: $output_file")
    end

    return fig
end
