using SurfaceWaterIntegratedModeling
export Domain3D, CellProperties, SimulationLayerSnapshot, ReservoirProperties

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
    shale_pressure_threshold::Float64
    brine_co2_density_difference::Float64
    residual_leakage_time::Float64

    function ReservoirProperties(sand_porosity::Float64 = 0.4,
                                 sand_residual_co2_saturation::Float64 = 0.2,
                                 sand_irreducible_water_saturation::Float64 = 0.3,
                                 shale_pressure_threshold::Float64 = 98000.0,
                                 brine_co2_density_difference::Float64 = 400.0,
                                 residual_leakage_time::Float64 = 1.0)
        new(sand_porosity,
            sand_residual_co2_saturation,
            sand_irreducible_water_saturation,
            shale_pressure_threshold,
            brine_co2_density_difference,
            residual_leakage_time)
    end
end

struct SimulationLayerSnapshot
    timestamp::Float64
    spill_event::SpillEvent
    heights::Array{Float64, 2}
    filled_traps::Vector{Bool}
    co2_volume::Float64
    mobile_co2_volume::Float64
    residual_trapped_co2_volume::Float64
    residual_trapped::Vector{Bool}
end

struct SimulationSnapshot
    timestamp::Float64
    total_co2_volume::Float64
    total_mobile_co2_volume::Float64
    total_residual_trapped_co2_volume::Float64
    layer_snapshots::Vector{SimulationLayerSnapshot}
end
