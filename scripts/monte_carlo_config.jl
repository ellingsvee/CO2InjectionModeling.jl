"""
Configuration for Monte Carlo uncertainty analysis.

Define parameter ranges and distributions for reservoir properties.
Each parameter can have a uniform, normal, or fixed distribution.
"""

using Distributions

struct ParameterDistribution
    name::String
    distribution::Distribution
end

"""
    MonteCarloConfig

Configuration for Monte Carlo simulation including:
- Parameter distributions for reservoir properties
- Simulation settings (number of realizations, time range, etc.)
"""
struct MonteCarloConfig
    # Number of Monte Carlo realizations
    n_realizations::Int

    # Parameter distributions
    param_distributions::Dict{Symbol, Distribution}

    # Simulation time settings
    start_time::Float64
    end_time::Float64
    num_snapshots::Int

    # Random seed for reproducibility
    random_seed::Union{Int, Nothing}
end

"""
    create_default_monte_carlo_config(; kwargs...)

Create a default Monte Carlo configuration with reasonable parameter ranges
based on Sleipner reservoir properties.

# Arguments
- `n_realizations::Int=100`: Number of Monte Carlo realizations
- `start_time::Float64=0.0`: Simulation start time (years)
- `end_time::Float64=15.0`: Simulation end time (years)
- `num_snapshots::Int=30`: Number of time snapshots
- `random_seed::Union{Int,Nothing}=nothing`: Random seed for reproducibility

# Returns
- `MonteCarloConfig`: Configuration object with parameter distributions
"""
function create_default_monte_carlo_config(;
    n_realizations::Int=100,
    start_time::Float64=0.0,
    end_time::Float64=15.0,
    num_snapshots::Int=30,
    random_seed::Union{Int,Nothing}=nothing
)
    # Define parameter distributions
    # Based on typical Sleipner reservoir properties with ±20% uncertainty

    param_distributions = Dict{Symbol, Distribution}(
        # Porosity: uniform distribution between 0.32 and 0.48 (0.4 ± 20%)
        :sand_porosity => Uniform(0.32, 0.48),

        # Residual CO2 saturation: uniform between 0.15 and 0.25 (0.2 ± 25%)
        :sand_residual_co2_saturation => Uniform(0.15, 0.25),

        # Irreducible water saturation: uniform between 0.24 and 0.36 (0.3 ± 20%)
        :sand_irreducible_water_saturation => Uniform(0.24, 0.36),

        # Shale pressure threshold: normal distribution around 98000 Pa with 10% std
        :shale_pressure_threshold => Normal(98000.0, 9800.0),

        # Leakage height: uniform between 16 and 24 m (20 ± 20%)
        :leakage_height => Uniform(16.0, 24.0),

        # Residual leakage time: uniform between 3 and 7 years (5 ± 40%)
        :residual_leakage_time => Uniform(3.0, 7.0)
    )

    return MonteCarloConfig(
        n_realizations,
        param_distributions,
        start_time,
        end_time,
        num_snapshots,
        random_seed
    )
end

"""
    create_custom_monte_carlo_config(; kwargs...)

Create a custom Monte Carlo configuration with user-specified parameter ranges.

# Arguments
- `n_realizations::Int=100`: Number of Monte Carlo realizations
- `sand_porosity::Tuple{Float64,Float64}=(0.32, 0.48)`: Min/max porosity
- `sand_residual_co2_saturation::Tuple{Float64,Float64}=(0.15, 0.25)`: Min/max residual saturation
- `sand_irreducible_water_saturation::Tuple{Float64,Float64}=(0.24, 0.36)`: Min/max water saturation
- `shale_pressure_threshold::Tuple{Float64,Float64}=(80000.0, 116000.0)`: Min/max pressure threshold
- `leakage_height::Tuple{Float64,Float64}=(16.0, 24.0)`: Min/max leakage height
- `residual_leakage_time::Tuple{Float64,Float64}=(3.0, 7.0)`: Min/max leakage time
- `start_time::Float64=0.0`: Simulation start time (years)
- `end_time::Float64=15.0`: Simulation end time (years)
- `num_snapshots::Int=30`: Number of time snapshots
- `random_seed::Union{Int,Nothing}=nothing`: Random seed

# Returns
- `MonteCarloConfig`: Configuration object
"""
function create_custom_monte_carlo_config(;
    n_realizations::Int=100,
    sand_porosity::Tuple{Float64,Float64}=(0.32, 0.48),
    sand_residual_co2_saturation::Tuple{Float64,Float64}=(0.15, 0.25),
    sand_irreducible_water_saturation::Tuple{Float64,Float64}=(0.24, 0.36),
    shale_pressure_threshold::Tuple{Float64,Float64}=(80000.0, 116000.0),
    leakage_height::Tuple{Float64,Float64}=(16.0, 24.0),
    residual_leakage_time::Tuple{Float64,Float64}=(3.0, 7.0),
    start_time::Float64=0.0,
    end_time::Float64=15.0,
    num_snapshots::Int=30,
    random_seed::Union{Int,Nothing}=nothing
)
    param_distributions = Dict{Symbol, Distribution}(
        :sand_porosity => Uniform(sand_porosity[1], sand_porosity[2]),
        :sand_residual_co2_saturation => Uniform(sand_residual_co2_saturation[1], sand_residual_co2_saturation[2]),
        :sand_irreducible_water_saturation => Uniform(sand_irreducible_water_saturation[1], sand_irreducible_water_saturation[2]),
        :shale_pressure_threshold => Uniform(shale_pressure_threshold[1], shale_pressure_threshold[2]),
        :leakage_height => Uniform(leakage_height[1], leakage_height[2]),
        :residual_leakage_time => Uniform(residual_leakage_time[1], residual_leakage_time[2])
    )

    return MonteCarloConfig(
        n_realizations,
        param_distributions,
        start_time,
        end_time,
        num_snapshots,
        random_seed
    )
end
