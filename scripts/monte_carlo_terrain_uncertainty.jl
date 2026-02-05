"""
Monte Carlo simulation with both reservoir parameter AND terrain uncertainty.

This script extends the standard Monte Carlo analysis by also varying the
topography/terrain for each realization. Since the trap structure depends on
topography, we must re-run spillanalysis (via analyze_base_surfaces) for each
perturbed topography.

Key differences from monte_carlo_simulation.jl:
1. Topography is perturbed for each realization
2. analyze_base_surfaces is called per-realization (computationally more expensive)
3. Injection events are regenerated for each new layer structure
"""

using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using CairoMakie
using Printf
using Random
using Distributions
using ProgressMeter

# Include helper modules
include("terrain_perturbation.jl")
include("monte_carlo_analysis.jl")


# =============================================================================
# TERRAIN UNCERTAINTY CONFIGURATION
# =============================================================================

"""
    TerrainUncertaintyConfig

Configuration for terrain/topography uncertainty in Monte Carlo simulation.
"""
struct TerrainUncertaintyConfig
    # Standard deviation of depth perturbation (meters)
    noise_std::Float64

    # Grid spacing for noise sampling (controls smoothness/wavelength)
    # Lower values = more high-frequency variation, higher = smoother perturbations
    sample_spacing::Int

    # Whether to perturb all surfaces or only sand layer tops
    perturb_all_surfaces::Bool

    # Correlation of noise between adjacent layers (0-1)
    # 0 = independent perturbations, 1 = identical perturbations
    correlation_between_layers::Float64
end

"""
    create_default_terrain_config(; kwargs...)

Create default terrain uncertainty configuration.

# Arguments
- `noise_std::Float64=1.0`: Standard deviation of depth perturbation (meters)
- `sample_spacing::Int=3`: Grid spacing for noise sampling
- `perturb_all_surfaces::Bool=true`: Perturb all surfaces
- `correlation_between_layers::Float64=0.5`: Vertical correlation of perturbations
"""
function create_default_terrain_config(;
    noise_std::Float64=1.0,
    sample_spacing::Int=3,
    perturb_all_surfaces::Bool=true,
    correlation_between_layers::Float64=0.5
)
    return TerrainUncertaintyConfig(
        noise_std,
        sample_spacing,
        perturb_all_surfaces,
        correlation_between_layers
    )
end


# =============================================================================
# COMBINED MONTE CARLO CONFIG
# =============================================================================

"""
    TerrainMonteCarloConfig

Combined configuration for Monte Carlo simulation with terrain uncertainty.
"""
struct TerrainMonteCarloConfig
    # Standard Monte Carlo config for reservoir parameters
    param_config::MonteCarloConfig

    # Terrain uncertainty config
    terrain_config::TerrainUncertaintyConfig
end


# =============================================================================
# TERRAIN-AWARE MONTE CARLO RUNNER
# =============================================================================

"""
    TerrainMonteCarloResult

Stores results from a single Monte Carlo realization with terrain uncertainty.
Extends MonteCarloResult with information about the perturbed topography.
"""
struct TerrainMonteCarloResult
    realization_id::Int
    parameters::Dict{Symbol, Float64}
    snapshots::Vector{ReservoirSnapshot}
    seqs::Union{Vector{Vector{SpillEvent}}, Nothing}
    leakage_states::Union{Vector{LeakageState}, Nothing}
    # Store some summary info about the perturbed topography (not the full arrays)
    topo_depth_range::Tuple{Float64, Float64}
end

"""
    TerrainMonteCarloEnsemble

Aggregated results from all Monte Carlo realizations with terrain uncertainty.
"""
struct TerrainMonteCarloEnsemble
    config::TerrainMonteCarloConfig
    results::Vector{TerrainMonteCarloResult}
    timepoints::Vector{Float64}
end


"""
    run_single_terrain_realization(
        realization_id::Int,
        config::TerrainMonteCarloConfig,
        base_topography::SleipnerTopography,
        boundary_condition::Symbol;
        store_seqs::Bool=false
    )

Run a single Monte Carlo realization with perturbed terrain.

# Arguments
- `realization_id::Int`: Identifier for this realization
- `config::TerrainMonteCarloConfig`: Combined configuration
- `base_topography::SleipnerTopography`: Original topography to perturb
- `boundary_condition::Symbol`: Boundary condition (:open or :closed)
- `store_seqs::Bool=false`: Whether to store spill sequences

# Returns
- `TerrainMonteCarloResult`: Results from this realization
"""
function run_single_terrain_realization(
    realization_id::Int,
    config::TerrainMonteCarloConfig,
    base_topography::SleipnerTopography,
    boundary_condition::Symbol;
    store_seqs::Bool=false
)
    tc = config.terrain_config
    pc = config.param_config

    # Step 1: Create perturbed topography
    perturbed_topography = create_perturbed_topography(
        base_topography;
        noise_std=tc.noise_std,
        sample_spacing=tc.sample_spacing,
        perturb_all_surfaces=tc.perturb_all_surfaces,
        correlation_between_layers=tc.correlation_between_layers
    )

    # Step 2: Create domain from perturbed topography
    domain = create_domain_from_topography(perturbed_topography, 1.0)

    # Step 3: Analyze trap structure for perturbed topography
    # This is the key difference - we recompute spillanalysis for each realization
    layers = analyze_base_surfaces(perturbed_topography; boundary_condition=boundary_condition)

    # Step 4: Get injection location (use same location relative to grid)
    # Note: The feeder location is absolute UTM, so grid index should be the same
    xy, _ = load_feeder_location(perturbed_topography)

    # Step 5: Generate injection events for the new layer structure
    injection_events = generate_sleipner_injection_events(layers, xy)

    # Step 6: Sample reservoir properties
    n_layers = length(layers)
    reservoir_properties, sampled_params = sample_reservoir_properties(pc, n_layers)

    # Step 7: Run simulation
    seqs, leakage_states = fill_layers(
        layers,
        domain,
        reservoir_properties,
        injection_events;
        verbose=false
    )

    # Step 8: Generate snapshots
    snapshots = generate_reservoir_snapshots(
        layers,
        seqs,
        leakage_states,
        domain,
        reservoir_properties,
        injection_events;
        num_snapshots=pc.num_snapshots,
        start_time=pc.start_time,
        end_time=pc.end_time,
        verbose=false
    )

    # Store optional data
    stored_seqs = store_seqs ? seqs : nothing
    stored_leakage_states = store_seqs ? leakage_states : nothing

    return TerrainMonteCarloResult(
        realization_id,
        sampled_params,
        snapshots,
        stored_seqs,
        stored_leakage_states,
        (perturbed_topography.depth_min, perturbed_topography.depth_max)
    )
end


"""
    run_terrain_monte_carlo_simulation(
        config::TerrainMonteCarloConfig,
        base_topography::SleipnerTopography;
        boundary_condition::Symbol=:closed,
        verbose::Bool=true,
        store_seqs::Bool=false
    )

Run Monte Carlo simulation with terrain uncertainty.

# Arguments
- `config::TerrainMonteCarloConfig`: Combined configuration
- `base_topography::SleipnerTopography`: Original topography to perturb
- `boundary_condition::Symbol=:closed`: Boundary condition for spillanalysis
- `verbose::Bool=true`: Show progress bar
- `store_seqs::Bool=false`: Store spill sequences for animations

# Returns
- `TerrainMonteCarloEnsemble`: Ensemble results
"""
function run_terrain_monte_carlo_simulation(
    config::TerrainMonteCarloConfig,
    base_topography::SleipnerTopography;
    boundary_condition::Symbol=:closed,
    verbose::Bool=true,
    store_seqs::Bool=false
)
    pc = config.param_config

    # Set random seed if specified
    if pc.random_seed !== nothing
        Random.seed!(pc.random_seed)
    end

    # Compute timepoints
    timepoints = range(pc.start_time, pc.end_time, length=pc.num_snapshots)

    # Initialize results storage
    results = Vector{TerrainMonteCarloResult}(undef, pc.n_realizations)

    # Progress bar
    progress = Progress(
        pc.n_realizations,
        desc="Running Monte Carlo with terrain uncertainty: ",
        enabled=verbose
    )

    # Run simulations
    for i in 1:pc.n_realizations
        results[i] = run_single_terrain_realization(
            i, config, base_topography, boundary_condition;
            store_seqs=store_seqs
        )
        next!(progress)
    end

    return TerrainMonteCarloEnsemble(config, results, collect(timepoints))
end


# =============================================================================
# CONVERSION FUNCTIONS FOR REUSING ANALYSIS CODE
# =============================================================================

"""
Convert TerrainMonteCarloEnsemble to MonteCarloEnsemble for analysis functions.
"""
function to_standard_ensemble(terrain_ensemble::TerrainMonteCarloEnsemble)
    # Convert TerrainMonteCarloResult to MonteCarloResult
    standard_results = [
        MonteCarloResult(
            r.realization_id,
            r.parameters,
            r.snapshots,
            r.seqs,
            r.leakage_states
        )
        for r in terrain_ensemble.results
    ]

    return MonteCarloEnsemble(
        terrain_ensemble.config.param_config,
        standard_results,
        terrain_ensemble.timepoints
    )
end


# =============================================================================
# MAIN SIMULATION SCRIPT
# =============================================================================


println("\n" * "="^80)
println("MONTE CARLO SIMULATION WITH TERRAIN UNCERTAINTY")
println("="^80)

# =========================================================================
# CONFIGURATION
# =========================================================================

# Parameter uncertainty configuration
param_config = create_default_monte_carlo_config(
n_realizations=25,
    start_time=0.0,
    end_time=15.0,
    num_snapshots=45,
    random_seed=42
)

# Terrain uncertainty configuration
terrain_config = create_default_terrain_config(
    noise_std=1.5,              # 3.0 meters standard deviation
    sample_spacing=3,           # Sample every 3rd grid point (~150m wavelength)
    perturb_all_surfaces=true,  # Perturb all depth surfaces
    correlation_between_layers=0.5  # 50% correlation between layers
)

# Combined configuration
config = TerrainMonteCarloConfig(param_config, terrain_config)

# =========================================================================
# SETUP
# =========================================================================

println("\nLoading base topography...")
boundary_condition = :closed
base_topography = load_sleipner_topography()

# =========================================================================
# RUN SIMULATION
# =========================================================================

terrain_ensemble = run_terrain_monte_carlo_simulation(
    config,
    base_topography;
    boundary_condition=boundary_condition,
    verbose=true,
    store_seqs=true
);

# Convert to standard ensemble for analysis
ensemble = to_standard_ensemble(terrain_ensemble);

# =========================================================================
# ANALYSIS AND VISUALIZATION
# =========================================================================

# Create output directory
output_dir = "monte_carlo_terrain_results"
mkpath(output_dir)

# Layer distribution barplot
plot_layer_distribution_barplot(
    ensemble;
    output_file="$output_dir/layer_distribution_barplot.svg",
    title=nothing,
    bar_color=colorant"#386624",
    bar_alpha=1.0,
    error_bar_color=:black,
    error_bar_linewidth=2.0,
    fontsize_labels=16,
    fontsize_ticks=15,
    fontsize_values=14,
    figure_size=(800, 300)
)

# Print summary statistics
print_percentages(ensemble)
