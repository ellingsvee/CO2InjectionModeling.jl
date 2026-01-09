using SurfaceWaterIntegratedModeling
export Domain3D, CellProperties, SimulationLayerSnapshot, SimulationSnapshot, ReservoirProperties
export LeakageRecord, LeakageState

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

    # Constructor with automatic calculation of cell sizes
    function Domain3D(nx, ny, nz, length_x, length_y, length_z, depth_min, depth_max)
        dx = length_x / nx
        dy = length_y / ny
        dz = length_z / nz
        new(nx, ny, nz, length_x, length_y, length_z, dx, dy, dz, depth_min, depth_max)
    end
end


struct InjectionEvent
    timestamp::Float64
    injection_rate::Union{Matrix{Float64}, Float64} 
end


struct CellProperties
    porosity::Array{Float64, 3}
    pressure_threshold::Array{Float64, 3}
    residual_co2_saturation::Array{Float64, 3}
    irreducible_water_saturation::Array{Float64, 3}

    function CellProperties(domain::Domain3D)
        nx, ny, nz = domain.nx, domain.ny, domain.nz
        new(zeros(nx, ny, nz), zeros(nx, ny, nz), zeros(nx, ny, nz), zeros(nx, ny, nz))
    end
end

struct ReservoirProperties
    sand_porosity::Float64
    sand_residual_co2_saturation::Float64
    sand_irreducible_water_saturation::Float64
    shale_pressure_threshold::Float64  # Mean pressure threshold (Pa)
    leakage_height::Float64  # Mean leakage height derived from pressure (m)
    residual_leakage_time::Float64
    # Spatial variability parameters (pressure-based)
    shale_pressure_threshold_std::Float64  # Std dev for spatial variation (Pa)
    shale_pressure_threshold_per_trap::Union{Vector{Float64}, Nothing}  # Pre-specified per-trap pressures
    brine_density::Float64  # kg/m³
    co2_density::Float64    # kg/m³
end

# Outer constructor with automatic leakage height computation
function ReservoirProperties(
    sand_porosity::Float64,
    sand_residual_co2_saturation::Float64,
    sand_irreducible_water_saturation::Float64,
    shale_pressure_threshold::Float64,
    residual_leakage_time::Float64;
    leakage_height::Union{Float64, Nothing}=nothing,
    shale_pressure_threshold_std::Float64=0.0,
    shale_pressure_threshold_per_trap::Union{Vector{Float64}, Nothing}=nothing,
    brine_density::Float64=1020.0,
    co2_density::Float64=460.0
)
    # Compute leakage height from pressure if not provided
    if isnothing(leakage_height)
        if isfinite(shale_pressure_threshold) && shale_pressure_threshold > 0.0
            g = 9.81  # m/s²
            density_diff = brine_density - co2_density
            leakage_height = shale_pressure_threshold / (density_diff * g)
        else
            leakage_height = Inf  # Impermeable caprock
        end
    end

    ReservoirProperties(
        sand_porosity,
        sand_residual_co2_saturation,
        sand_irreducible_water_saturation,
        shale_pressure_threshold,
        leakage_height,
        residual_leakage_time,
        shale_pressure_threshold_std,
        shale_pressure_threshold_per_trap,
        brine_density,
        co2_density
    )
end

# Backward compatibility: old constructor with explicit leakage_height
function ReservoirProperties(
    sand_porosity::Float64,
    sand_residual_co2_saturation::Float64,
    sand_irreducible_water_saturation::Float64,
    shale_pressure_threshold::Float64,
    leakage_height::Float64,
    residual_leakage_time::Float64;
    shale_pressure_threshold_std::Float64=0.0,
    shale_pressure_threshold_per_trap::Union{Vector{Float64}, Nothing}=nothing,
    brine_density::Float64=1020.0,
    co2_density::Float64=460.0
)
    ReservoirProperties(
        sand_porosity,
        sand_residual_co2_saturation,
        sand_irreducible_water_saturation,
        shale_pressure_threshold,
        residual_leakage_time;
        leakage_height=leakage_height,
        shale_pressure_threshold_std=shale_pressure_threshold_std,
        shale_pressure_threshold_per_trap=shale_pressure_threshold_per_trap,
        brine_density=brine_density,
        co2_density=co2_density
    )
end

struct SimulationLayerSnapshot
    timestamp::Float64
    spill_event::SpillEvent
    filled_traps::Vector{Bool}
    injected_volume::Float64
    co2_volume::Float64
end

"""
    SimulationSnapshot

Represents the state of the entire simulation at a point in time.
Contains snapshots from all layers.
"""
struct SimulationSnapshot
    timestamp::Float64
    total_injected_volume::Float64
    total_co2_volume::Float64
    layer_snapshots::Vector{SimulationLayerSnapshot}
end

"""
    LeakageRecord

Records a leakage event for generating WeatherEvents in the overlying layer.
- `start_time`: When leakage started
- `trap_id`: The trap that is leaking
- `leakage_location`: Grid cell (CartesianIndex) where leakage occurs (trap's lowest point)
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
- `residual_saturation`: Fraction of CO2 that remains after drainage (from ReservoirProperties)
- `residual_leakage_time`: Time over which residual drainage occurs (from ReservoirProperties)
"""
mutable struct LeakageState
    leaking::Vector{Bool}
    draining::Vector{Bool}
    leakage_volume::Vector{Float64}
    leakage_start_time::Vector{Float64}
    leakage_records::Vector{LeakageRecord}
    leakage_height::Vector{Float64}  # Per-trap heights
    # Residual leakage fields
    initial_volume_at_leak::Vector{Float64}
    residual_saturation::Float64
    residual_leakage_time::Float64
end
