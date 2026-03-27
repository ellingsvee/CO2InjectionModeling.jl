using SurfaceWaterIntegratedModeling
export Domain3D, CellProperties, SimulationLayerSnapshot, SimulationSnapshot, ReservoirProperties
export LeakageRecord, LeakageState
export InjectionEvent
export LeakageState, LeakageRecord

"""
    Domain3D

3D domain describing the reservoir grid.

# Fields
- `nx`, `ny`, `nz`: Grid cell counts in x, y, and z directions
- `length_x`, `length_y`, `length_z`: Physical domain extent in meters
- `dx`, `dy`, `dz`: Cell sizes in meters (computed from length / n)
- `depth_min`, `depth_max`: Shallowest and deepest depths in the domain (m)

Construct via [`create_domain`](@ref) from a topography object, or directly:

    Domain3D(nx, ny, nz, length_x, length_y, length_z, depth_min, depth_max)
"""
struct Domain3D
    nx::Int
    ny::Int
    nz::Int
    length_x::Float64
    length_y::Float64
    length_z::Float64
    dx::Float64
    dy::Float64
    dz::Float64
    depth_min::Float64
    depth_max::Float64

    function Domain3D(nx, ny, nz, length_x, length_y, length_z, depth_min, depth_max)
        dx = length_x / nx
        dy = length_y / ny
        dz = length_z / nz
        new(nx, ny, nz, length_x, length_y, length_z, dx, dy, dz, depth_min, depth_max)
    end
end


"""
    InjectionEvent(timestamp, injection_rate)

Defines a CO2 injection rate starting at `timestamp`.

The simulation treats injection as piecewise-constant: a rate set at `timestamp`
remains in effect until the next `InjectionEvent`.  To stop injection, add an
event with `injection_rate` equal to zero.

# Fields
- `timestamp`: Time (in model time units) at which this rate becomes active
- `injection_rate`: Physical injection rate in m^3/year per grid cell,
  either a scalar or a matrix of size `(nx + 2*pad, ny + 2*pad)`

# Example
```julia
pad = 2
rate = zeros(NX + 2pad, NY + 2pad)
rate[div(NX,2)+pad, div(NY,2)+pad] = 25_000.0  # m^3/yr at centre cell
events = [InjectionEvent(0.0, rate), InjectionEvent(10.0, zeros(size(rate)))]
```
"""
struct InjectionEvent
    timestamp::Float64
    injection_rate::Union{Matrix{Float64},Float64}
end


struct CellProperties
    porosity::Array{Float64,3}
    pressure_threshold::Array{Float64,3}
    residual_co2_saturation::Array{Float64,3}
    irreducible_water_saturation::Array{Float64,3}

    function CellProperties(domain::Domain3D)
        nx, ny, nz = domain.nx, domain.ny, domain.nz
        new(zeros(nx, ny, nz), zeros(nx, ny, nz), zeros(nx, ny, nz), zeros(nx, ny, nz))
    end
end

"""
    ReservoirProperties(sand_porosity, sand_residual_co2_saturation,
                        sand_irreducible_water_saturation,
                        shale_pressure_threshold, residual_leakage_time;
                        brine_density=1020.0, co2_density=460.0)

Physical properties of the reservoir rock and fluids.

# Arguments
- `sand_porosity`: Porosity of the sand (0–1)
- `sand_residual_co2_saturation`: Residual (irreducible) CO2 saturation after drainage (0–1)
- `sand_irreducible_water_saturation`: Irreducible water saturation in sand (0–1)
- `shale_pressure_threshold`: Caprock entry pressure in Pa (`Inf` for impermeable caprock)
- `residual_leakage_time`: Duration over which residual CO2 drains after leakage onset (years)
- `brine_density`: Brine density in kg/m^3 (default 1020.0)
- `co2_density`: CO2 density in kg/m^3 (default 460.0)

# Derived fields
- `leakage_height`: CO2 column height threshold that triggers leakage (m),
  computed as `shale_pressure_threshold / ((brine_density - co2_density) * g)`.
  Set to `Inf` when the caprock is impermeable.

# Example
```julia
rp = ReservoirProperties(0.3, 0.4, 0.1, 15_000.0, 5.0)
rp_cap = ReservoirProperties(0.3, 0.4, 0.1, Inf, 5.0)  # sealed caprock
```
"""
struct ReservoirProperties
    sand_porosity::Float64
    sand_residual_co2_saturation::Float64
    sand_irreducible_water_saturation::Float64
    shale_pressure_threshold::Union{Float64,Vector{Float64}}  # Mean pressure threshold (Pa)
    leakage_height::Union{Float64,Vector{Float64}}  # Mean leakage height derived from pressure (m)
    residual_leakage_time::Float64
    brine_density::Float64  # kg/m^3
    co2_density::Float64    # kg/m^3

    function ReservoirProperties(
        sand_porosity::Float64,
        sand_residual_co2_saturation::Float64,
        sand_irreducible_water_saturation::Float64,
        shale_pressure_threshold::Union{Float64,Vector{Float64}},
        residual_leakage_time::Float64;
        brine_density::Float64=1020.0,
        co2_density::Float64=460.0
    )

        # Compute leakage height from pressure if not provided
        # pressure == 0.0: shale-breaks (leakage_height = 0.0, immediate pass-through)
        # pressure > 0.0: normal leakage threshold
        # pressure == Inf: impermeable caprock (leakage_height = Inf)
        g = 9.81
        density_diff = brine_density - co2_density
        if isa(shale_pressure_threshold, Float64)
            if shale_pressure_threshold == 0.0
                leakage_height = 0.0  # Shale-break: immediate pass-through
            elseif isfinite(shale_pressure_threshold) && shale_pressure_threshold > 0.0
                leakage_height = shale_pressure_threshold / (density_diff * g)
            else
                leakage_height = Inf  # Impermeable caprock
            end
        else
            leakage_height = similar(shale_pressure_threshold)
            for i in eachindex(shale_pressure_threshold)
                if shale_pressure_threshold[i] == 0.0
                    leakage_height[i] = 0.0  # Shale-break: immediate pass-through
                elseif isfinite(shale_pressure_threshold[i]) && shale_pressure_threshold[i] > 0.0
                    leakage_height[i] = shale_pressure_threshold[i] / (density_diff * g)
                else
                    leakage_height[i] = Inf  # Impermeable caprock
                end
            end
        end

        new(
            sand_porosity,
            sand_residual_co2_saturation,
            sand_irreducible_water_saturation,
            shale_pressure_threshold,
            leakage_height,
            residual_leakage_time,
            brine_density,
            co2_density
        )
    end

end

"""
LeakageRecord. Records a leakage event for generating WeatherEvents in the overlying layer.
"""
struct LeakageRecord
    start_time::Float64
    trap_id::Int
    leakage_location::CartesianIndex{2}
end

"""
    LeakageState

Tracks the leakage state for all traps in a layer during simulation.
- `leaking`: Boolean vector indicating if each trap has reached leakage threshold (edge=0)
- `draining`: Boolean vector indicating if each trap is experiencing residual drainage
  (includes leaking traps AND their filled ancestors whose CO2 drains through)
- `leakage_volume`: Volume at which leakage starts for each trap (precomputed from leakage_height)
- `leakage_start_time`: When leakage/drainage started for each trap (Inf if not yet)
- `leakage_records`: Vector of LeakageRecord for generating upstream WeatherEvents
- `leakage_height`: Per-trap threshold heights for leakage (sampled from ReservoirProperties)
- `initial_volume_at_leak`: Volume in each trap when drainage started (for residual drainage)
- `residual_volume_fraction`: Fraction of SWIM volume that remains after drainage, equal to
  `S_r / (1 - S_wi)` where `S_r` is residual CO2 saturation and `S_wi` is irreducible water
  saturation. This converts from pore-space saturation to the equivalent SWIM volume fraction.
- `residual_leakage_time`: Time over which residual drainage occurs (from ReservoirProperties)
- `cumulative_no_inflow_time`: Total time each leaking trap has spent without inflow (for dynamic equilibrium)
- `volume_at_last_state_change`: Stored volume when inflow state last changed
- `time_of_last_state_change`: Time when inflow state last changed
- `has_inflow`: Whether each leaking trap currently has positive inflow
"""
mutable struct LeakageState
    leaking::Vector{Bool}
    draining::Vector{Bool}
    leakage_volume::Vector{Float64}
    leakage_start_time::Vector{Float64}
    leakage_records::Vector{LeakageRecord}
    leakage_height::Vector{Float64}
    # Residual leakage fields
    initial_volume_at_leak::Vector{Float64}
    residual_volume_fraction::Float64
    residual_leakage_time::Float64
    # Dynamic equilibrium tracking (for directly leaking traps)
    cumulative_no_inflow_time::Vector{Float64}
    volume_at_last_state_change::Vector{Float64}
    time_of_last_state_change::Vector{Float64}
    has_inflow::Vector{Bool}
end
