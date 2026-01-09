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
        # :sand_residual_co2_saturation => Uniform(0.15, 0.25),
        :sand_residual_co2_saturation => Normal(0.2, 0),

        # Irreducible water saturation: uniform between 0.24 and 0.36 (0.3 ± 20%)
        # :sand_irreducible_water_saturation => Uniform(0.25, 0.35),
        :sand_irreducible_water_saturation => Normal(0.3, 0),

        # Shale pressure threshold: normal distribution around 98000 Pa with 10% std
        :shale_pressure_threshold => Normal(98000.0, 0),
        # :shale_pressure_threshold => Uniform(49000.0, 147000.0),

        # Residual leakage time: uniform between 3 and 7 years (5 ± 40%)
        :residual_leakage_time => Uniform(3.0, 7.0),

        # Shale pressure threshold spatial variability: std deviation in Pa (0 = no variability)
        # User specified 16000 Pa as the std for spatial variation
        # :shale_pressure_threshold_std => Normal(16000.0, 0)
        # :shale_pressure_threshold_std => Normal(8000.0, 0)
        :shale_pressure_threshold_std => Normal(0.0, 0)
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


