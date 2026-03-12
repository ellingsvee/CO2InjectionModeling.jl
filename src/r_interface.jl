"""
    Generic R interface for CO2BatchFill.

Provides dataset-agnostic functions for running CO2 injection simulations
from R (via JuliaCall) or any other external interface.

Key entry points:
- `setup_simulator(topography; ...)` — initialize from any AbstractTopography
- `setup_simulator_from_surfaces(; ...)` — initialize from raw surface arrays
- `configure_reservoir(; ...)` — set reservoir properties
- `run_simulation(; ...)` — run simulation and return results
- `generate_birdseye_animation(; ...)` / `generate_cross_section_animation(; ...)` — create animations
"""

using SurfaceWaterIntegratedModeling: SpillEvent, WeatherEvent, numtraps

export SimulatorState
export setup_simulator, setup_simulator_from_surfaces
export configure_reservoir, run_simulation
export generate_cross_section_animation, generate_birdseye_animation

# Global state to hold simulator configuration
mutable struct SimulatorState
    topography::Union{Nothing,AbstractTopography}
    domain::Union{Nothing,Domain3D}
    layers::Union{Nothing,Vector{Layer}}
    reservoir_properties::Union{Nothing,Vector{ReservoirProperties}}
    boundary_condition::Symbol
    last_snapshots::Union{Nothing,Vector{MultiLayerSnapshot}}
    last_seqs::Union{Nothing,Vector{Vector{SpillEvent}}}
    last_weather_events_per_layer::Union{Nothing,Vector{Vector{WeatherEvent}}}
end

const SIMULATOR = SimulatorState(nothing, nothing, nothing, nothing, :open, nothing, nothing, nothing)

"""
    setup_simulator(topography::AbstractTopography; boundary_condition="open")

Set up the simulator from any `AbstractTopography` implementation.

# Arguments
- `topography`: Any `AbstractTopography` (e.g., `GenericTopography`)
- `boundary_condition`: `"open"` or `"closed"` (default: `"open"`)

# Returns
Dictionary with `status`, `n_layers`, `nx`, `ny`, `boundary_condition`, `nx_after_bc`, `ny_after_bc`
"""
function setup_simulator(
    topography::AbstractTopography;
    boundary_condition::String="open",
)
    try
        bc_symbol = Symbol(boundary_condition)
        if !(bc_symbol in [:open, :closed])
            return Dict("status" => "error", "message" => "boundary_condition must be 'open' or 'closed'")
        end

        SIMULATOR.topography = topography
        SIMULATOR.domain = create_domain(topography, 1.0)
        SIMULATOR.layers = analyze_base_surfaces(topography; boundary_condition=bc_symbol)
        SIMULATOR.boundary_condition = bc_symbol
        SIMULATOR.reservoir_properties = nothing

        nx, ny = get_grid_size(SIMULATOR.layers)
        nx_after_bc, ny_after_bc = get_padded_grid_size(SIMULATOR.layers)

        return Dict(
            "status" => "success",
            "n_layers" => length(SIMULATOR.layers),
            "nx" => nx, "ny" => ny,
            "boundary_condition" => boundary_condition,
            "nx_after_bc" => nx_after_bc, "ny_after_bc" => ny_after_bc,
        )
    catch e
        return Dict("status" => "error", "message" => string(e))
    end
end

"""
    setup_simulator_from_surfaces(; layer_tops, layer_bases, layer_names, dx, dy,
                                    boundary_condition="open", caprock_surface=nothing)

Set up the simulator from raw surface arrays (primary entry point for R users with custom data).

# Arguments
- `layer_tops`: Vector of 2D matrices (each nx × ny), top surfaces for each layer
- `layer_bases`: Vector of 2D matrices (each nx × ny), base surfaces for each layer
- `layer_names`: Vector of layer name strings (e.g., [\"L1\", \"L2\", ...])
- `dx`: Grid spacing in x direction (meters)
- `dy`: Grid spacing in y direction (meters)
- `boundary_condition`: `"open"` or `"closed"` (default: `"open"`)
- `caprock_surface`: Optional caprock top surface matrix (nx × ny)

# Returns
Dictionary with `status`, `n_layers`, `nx`, `ny`, `boundary_condition`, `nx_after_bc`, `ny_after_bc`
"""
function setup_simulator_from_surfaces(;
    layer_tops::Vector,
    layer_bases::Vector,
    layer_names::Vector{String},
    dx::Float64,
    dy::Float64,
    boundary_condition::String="open",
    caprock_surface::Union{Matrix{Float64},Nothing}=nothing
)
    try
        n = length(layer_tops)
        if length(layer_bases) != n || length(layer_names) != n
            return Dict("status" => "error", "message" => "layer_tops, layer_bases, and layer_names must have the same length")
        end
        if n == 0
            return Dict("status" => "error", "message" => "Must provide at least one layer")
        end

        # Build sand layer dicts
        sand_layers = Dict{String,Any}[]
        for i in 1:n
            top = Matrix{Float64}(layer_tops[i])
            base = Matrix{Float64}(layer_bases[i])
            push!(sand_layers, Dict{String,Any}("name" => layer_names[i], "top" => top, "base" => base))
        end

        nx, ny = size(Matrix{Float64}(layer_tops[1]))
        depth_min = minimum(minimum(Matrix{Float64}(t)) for t in layer_tops)
        depth_max = maximum(maximum(Matrix{Float64}(b)) for b in layer_bases)

        topo = GenericTopography(sand_layers, nx, ny, dx, dy, depth_min, depth_max;
            caprock_surface=caprock_surface)

        return setup_simulator(topo; boundary_condition=boundary_condition)
    catch e
        return Dict("status" => "error", "message" => string(e))
    end
end

"""
    configure_reservoir(; porosity, residual_co2_sat, irreducible_water_sat,
                         shale_pressure_threshold, residual_leakage_time,
                         brine_density=1020.0, co2_density=460.0, layer_specific=false)

Configure reservoir properties for the simulation.

# Arguments
- `porosity`: Sand porosity (0-1). Scalar or vector of length n_layers
- `residual_co2_sat`: Residual CO2 saturation (0-1). Scalar or vector
- `irreducible_water_sat`: Irreducible water saturation (0-1). Scalar or vector
- `shale_pressure_threshold`: Shale pressure threshold (Pa). Scalar or vector. Use `Inf` for impermeable caprock
- `residual_leakage_time`: Residual leakage time (years). Scalar or vector
- `brine_density`: Brine density (kg/m³). Scalar or vector (default: 1020.0)
- `co2_density`: CO2 density (kg/m³). Scalar or vector (default: 460.0)
- `layer_specific`: Set to `true` to provide vectors for layer-specific properties (default: `false`)
"""
function configure_reservoir(;
    porosity::Union{Float64,Vector{Float64}},
    residual_co2_sat::Union{Float64,Vector{Float64}},
    irreducible_water_sat::Union{Float64,Vector{Float64}},
    shale_pressure_threshold::Union{Float64,Vector{Float64}},
    residual_leakage_time::Union{Float64,Vector{Float64}},
    brine_density::Union{Float64,Vector{Float64}}=1020.0,
    co2_density::Union{Float64,Vector{Float64}}=460.0,
    layer_specific::Bool=false
)
    try
        if isnothing(SIMULATOR.layers)
            return Dict("status" => "error", "message" => "Must call setup_simulator first")
        end

        n_layers = length(SIMULATOR.layers)

        # Expand any scalar parameters to vectors of length n_layers
        _expand(x, n) = x isa Float64 ? fill(x, n) : x
        porosity = _expand(porosity, n_layers)
        residual_co2_sat = _expand(residual_co2_sat, n_layers)
        irreducible_water_sat = _expand(irreducible_water_sat, n_layers)
        shale_pressure_threshold = _expand(shale_pressure_threshold, n_layers)
        residual_leakage_time = _expand(residual_leakage_time, n_layers)
        brine_density = _expand(brine_density, n_layers)
        co2_density = _expand(co2_density, n_layers)

        SIMULATOR.reservoir_properties = [
            ReservoirProperties(
                porosity[i], residual_co2_sat[i], irreducible_water_sat[i],
                shale_pressure_threshold[i], residual_leakage_time[i];
                brine_density=brine_density[i], co2_density=co2_density[i]
            )
            for i in 1:n_layers
        ]

        return Dict("status" => "success", "n_layers" => n_layers)
    catch e
        return Dict("status" => "error", "message" => string(e))
    end
end

"""
    run_simulation(; start_time, end_time, time_step, injection_rate_matrices, verbose=false)

Run a CO2 injection simulation with specified parameters.

# Arguments
- `start_time`: Simulation start time (years)
- `end_time`: Simulation end time (years)
- `time_step`: Time step for output snapshots (years)
- `injection_rate_matrices`: Vector of 3D arrays (one per layer), each `n_times × nx × ny`
- `verbose`: Print progress messages (default: `false`)
"""
function run_simulation(;
    start_time::Float64,
    end_time::Float64,
    time_step::Float64,
    injection_rate_matrices::Vector,
    verbose::Bool=false
)
    try
        if isnothing(SIMULATOR.layers) || isnothing(SIMULATOR.reservoir_properties)
            return Dict("status" => "error", "message" => "Must call setup_simulator and configure_reservoir first")
        end

        n_layers = length(SIMULATOR.layers)
        grid_size = size(SIMULATOR.layers[1].trap_structure.topography)

        if length(injection_rate_matrices) != n_layers
            return Dict("status" => "error", "message" => "Must provide injection_rate_matrices with length $n_layers")
        end

        injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)

        for layer_idx in 1:n_layers
            rate_matrix = injection_rate_matrices[layer_idx]
            if ndims(rate_matrix) != 3
                return Dict("status" => "error", "message" => "Layer $layer_idx: injection_rate_matrix must be 3D")
            end

            n_times = size(rate_matrix, 1)
            layer_events = InjectionEvent[]
            for t_idx in 1:n_times
                time = start_time + (t_idx - 1) * time_step
                push!(layer_events, InjectionEvent(time, rate_matrix[t_idx, :, :]))
            end
            injection_events[layer_idx] = isempty(layer_events) ? [InjectionEvent(0.0, zeros(grid_size))] : layer_events
        end

        verbose && println("Running simulation from $start_time to $end_time years")

        seqs, leakage_states, weather_events_per_layer = fill_layers(
            SIMULATOR.layers, SIMULATOR.domain, SIMULATOR.reservoir_properties,
            injection_events; verbose=verbose
        )

        timepoints = collect(range(start_time, stop=end_time, step=time_step))

        snapshots = generate_multi_layer_snapshots(
            SIMULATOR.layers, seqs, leakage_states, weather_events_per_layer, timepoints
        )

        SIMULATOR.last_snapshots = snapshots
        SIMULATOR.last_seqs = seqs
        SIMULATOR.last_weather_events_per_layer = weather_events_per_layer

        rp = SIMULATOR.reservoir_properties
        domain = SIMULATOR.domain

        total_co2_volumes = [
            swim_volume_to_physical_volume(s.total_stored, rp[1], domain)
            for s in snapshots
        ]
        layer_co2_volumes = hcat([
            [swim_volume_to_physical_volume(s.layers[i].total_stored, rp[i], domain) for i in 1:n_layers]
            for s in snapshots
        ]...)'

        return Dict(
            "status" => "success",
            "timepoints" => [s.timestamp for s in snapshots],
            "total_co2_volumes" => total_co2_volumes,
            "layer_co2_volumes" => layer_co2_volumes,
            "num_layers" => n_layers,
            "num_traps_per_layer" => [numtraps(l.trap_structure) for l in SIMULATOR.layers]
        )
    catch e
        return Dict("status" => "error", "message" => string(e), "stacktrace" => sprint(showerror, e, catch_backtrace()))
    end
end

"""
    generate_cross_section_animation(; kwargs...)

Generate animation of CO2 trap filling from the last simulation.
"""
function generate_cross_section_animation(;
    output_file::String="multi_layer_filling.gif",
    num_frames::Int=30,
    start_time::Union{Float64,Nothing}=nothing,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=2,
    colormap::String="thermal",
    max_co2_height::Float64=20.0,
    show_contours::Bool=true,
    contour_levels::Int=10
)
    try
        if isnothing(SIMULATOR.last_seqs)
            return Dict("status" => "error", "message" => "Must call run_simulation first")
        end
        animate_multi_layer_filling(
            SIMULATOR.layers, SIMULATOR.last_seqs, SIMULATOR.domain;
            output_file=output_file, num_frames=num_frames, start_time=start_time,
            end_time=end_time, fps=fps, colormap=Symbol(colormap), max_co2_height=max_co2_height,
            show_contours=show_contours, contour_levels=contour_levels
        )
        return Dict("status" => "success", "output_file" => output_file)
    catch e
        return Dict("status" => "error", "message" => string(e))
    end
end

"""
    generate_birdseye_animation(; kwargs...)

Generate bird's eye view animation of CO2 trap filling.
"""
function generate_birdseye_animation(;
    output_file::String="multi_layer_filling.gif",
    num_frames::Int=30,
    start_time::Union{Float64,Nothing}=nothing,
    end_time::Union{Float64,Nothing}=nothing,
    fps::Int=2,
    colormap::String="thermal",
    max_co2_height::Float64=20.0,
    show_contours::Bool=true,
    contour_levels::Int=10
)
    generate_cross_section_animation(;
        output_file=output_file, num_frames=num_frames, start_time=start_time,
        end_time=end_time, fps=fps, colormap=colormap, max_co2_height=max_co2_height,
        show_contours=show_contours, contour_levels=contour_levels
    )
end
