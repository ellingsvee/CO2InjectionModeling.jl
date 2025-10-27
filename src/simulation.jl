# Simulation utilities for CO2 injection modeling

using SurfaceWaterIntegratedModeling
using Statistics

export compute_leakage_heights,
    interpolate_simulation_results,
    run_ensemble_simulations,
    compute_ensemble_statistics

"""
    compute_leakage_heights(tstructs; fraction=1.5)

Compute leakage heights for each layer based on spillpoint elevations.

The leakage height is set to `max_spillpoint_height / fraction` for each layer,
where `max_spillpoint_height` is the maximum height of spillpoints above trap bottoms.

# Arguments
- `tstructs::Vector{TrapStructure}`: Trap structures for all layers
- `fraction::Float64=1.5`: Divisor for max spillpoint height (higher = lower leakage threshold)

# Returns
- `Vector{Float64}`: Leakage heights for each layer (length = n_layers - 1)
"""
function compute_leakage_heights(tstructs::Vector{<:TrapStructure}; fraction::Float64=1.5)
    leakage_heights = Float64[]

    for i in 1:(length(tstructs)-1)
        tstruct = tstructs[i]
        spillpoints = tstruct.spillpoints
        max_leakage_height = -Inf

        for (j, spillpoint) in enumerate(spillpoints)
            trap_bottom = minimum(tstruct.topography[tstruct.footprints[j]])
            leakage_height = spillpoint.elevation - trap_bottom
            max_leakage_height = max(max_leakage_height, leakage_height)
        end

        push!(leakage_heights, max_leakage_height / fraction)
    end

    return leakage_heights
end

"""
    interpolate_simulation_results(tstructs, seqs, times_at_leakage, times; 
                                   filled_color=2, trap_color=8, river_color=6)

Interpolate CO2 distribution at specified time points for all layers.

# Arguments
- `tstructs::Vector{TrapStructure}`: Trap structures for all layers
- `seqs`: Simulation sequences for all layers
- `times_at_leakage::Vector{Float64}`: Times when CO2 leaks from each layer (filtered, no nothing values)
- `times::Vector{Float64}`: Time points at which to interpolate
- `filled_color::Int=2`: Color code for filled state
- `trap_color::Int=8`: Color code for trapped state
- `river_color::Int=6`: Color code for river state

# Returns
- `Vector{Vector}`: Interpolated textures for each layer at each time point
- `Int`: Number of layers that actually contain CO2
"""
function interpolate_simulation_results(
    tstructs::Vector{<:TrapStructure},
    seqs,
    times_at_leakage::Vector{Float64},
    times::Vector{Float64};
    filled_color::Int=2,
    trap_color::Int=1,
    river_color::Int=2
)
    # Determine how many layers actually have CO2
    n_filled_layers = length(times_at_leakage) + 1

    # Interpolate CO2 distribution for each layer
    texs_for_times = []

    for i in 1:length(tstructs)
        tstruct = tstructs[i]
        seq = seqs[i]

        # Check if this layer has any CO2
        if i > n_filled_layers
            # This layer never receives CO2 - create empty texture for all times
            empty_state = ones(Int, size(tstruct.topography))  # 1 = empty state
            tex = [empty_state for _ in times]
            push!(texs_for_times, tex)
            continue
        end

        seq_start = seq[1].timestamp
        seq_end = seq[end].timestamp

        # Determine active time range for this layer
        if i <= length(times_at_leakage)
            # This layer has leakage - active until leakage time
            layer_end_time = times_at_leakage[i]
        else
            # This is the last filled layer - active through sequence end
            layer_end_time = seq_end
        end

        # Get all requested times, clamped to sequence bounds
        times_bounded = [clamp(t, seq_start, seq_end) for t in times]
        times_for_interp = unique(sort(times_bounded))

        # Ensure we have at least 2 points for interpolation
        if length(times_for_interp) < 2
            times_for_interp = [seq_start, seq_end]
        end

        # Interpolate for all requested times
        tex_interp, = interpolate_timeseries(tstruct, seq, times_for_interp;
            filled_color=filled_color,
            trap_color=trap_color,
            river_color=river_color)

        # Map interpolated results back to original time grid
        tex = []
        for t in times
            t_clamped = clamp(t, seq_start, seq_end)
            idx = findfirst(x -> x == t_clamped, times_for_interp)
            if idx !== nothing
                push!(tex, tex_interp[idx])
            else
                # Shouldn't happen, but use final state as fallback
                push!(tex, tex_interp[end])
            end
        end

        push!(texs_for_times, tex)
    end

    return texs_for_times, n_filled_layers
end

"""
    run_ensemble_simulations(tstructs, injection_location, injection_rate, 
                            leakage_height_samples; n_common_times=100, start_time=0.1)

Run multiple simulations with different leakage heights and interpolate results 
onto a common time grid for comparison.

# Arguments
- `tstructs::Vector{TrapStructure}`: Trap structures for all layers
- `injection_location::Int`: Injection location in top layer
- `injection_rate::Float64`: Injection rate
- `leakage_height_samples::Vector{Vector{Float64}}`: Vector of leakage height vectors for each simulation
- `n_common_times::Int=100`: Number of time points for comparison
- `start_time::Float64=0.1`: Start time for common time grid

# Returns
- `common_times::Vector{Float64}`: Common time points for all simulations
- `all_texs::Vector{Vector{Vector}}`: Textures for each simulation, layer, and time
- `all_n_filled::Vector{Int}`: Number of filled layers for each simulation
- `simulation_metadata::Vector{NamedTuple}`: Metadata for each simulation (times_at_leakage, etc.)
"""
function run_ensemble_simulations(
    tstructs::Vector{<:TrapStructure},
    injection_location::Int,
    injection_rate::Float64,
    leakage_height_samples::Vector{Vector{Float64}};
    n_common_times::Int=100,
    start_time::Float64=0.1
)
    n_simulations = length(leakage_height_samples)

    # Storage for all simulation results
    all_seqs = []
    all_times_at_leakage = []
    all_leakage_locations = []

    println("Running $n_simulations ensemble simulations...")

    # Run all simulations
    for (sim_idx, leakage_heights) in enumerate(leakage_height_samples)
        println("\nSimulation $sim_idx/$n_simulations")
        seqs, times_at_leakage, leakage_locations = run_injection(
            tstructs, injection_location, injection_rate, leakage_heights
        )

        # Filter out nothing values
        times_at_leakage_filtered = Float64[x for x in times_at_leakage if x !== nothing]

        push!(all_seqs, seqs)
        push!(all_times_at_leakage, times_at_leakage_filtered)
        push!(all_leakage_locations, leakage_locations)
    end

    # Find maximum end time across all simulations
    max_end_time = -Inf
    for (sim_idx, seqs) in enumerate(all_seqs)
        n_filled = length(all_times_at_leakage[sim_idx]) + 1
        end_time = seqs[n_filled][end].timestamp
        max_end_time = max(max_end_time, end_time)
    end

    println("\nCommon time range: $start_time to $max_end_time")

    # Create common time grid
    common_times = collect(range(start_time, max_end_time, length=n_common_times))

    # Interpolate all simulations onto common time grid
    all_texs = []
    all_n_filled = Int[]

    println("\nInterpolating results onto common time grid...")
    for (sim_idx, seqs) in enumerate(all_seqs)
        texs, n_filled = interpolate_simulation_results(
            tstructs, seqs, all_times_at_leakage[sim_idx], common_times
        )
        push!(all_texs, texs)
        push!(all_n_filled, n_filled)
        println("Simulation $sim_idx: CO2 reached $n_filled layers")
    end

    # Create metadata for each simulation
    simulation_metadata = [
        (
            times_at_leakage=all_times_at_leakage[i],
            leakage_locations=all_leakage_locations[i],
            n_filled_layers=all_n_filled[i],
            leakage_heights=leakage_height_samples[i]
        )
        for i in 1:n_simulations
    ]

    return common_times, all_texs, all_n_filled, simulation_metadata
end

"""
    compute_ensemble_statistics(all_texs, layer_idx, time_idx)

Compute fraction of simulations with CO2 present at each location.

# Arguments
- `all_texs::Vector{Vector{Vector}}`: Textures from `run_ensemble_simulations`
- `layer_idx::Int`: Layer index
- `time_idx::Int`: Time index

# Returns
- `fraction_filled::Matrix{Float64}`: Fraction of simulations with CO2 present (0.0 to 1.0)
- `n_filled::Matrix{Int}`: Number of simulations with CO2 present at each location
"""
function compute_ensemble_statistics(all_texs::Vector, layer_idx::Int, time_idx::Int)
    n_simulations = length(all_texs)

    # Find a simulation that has data for this layer and time to get dimensions
    ny, nx = 0, 0
    for texs in all_texs
        if layer_idx <= length(texs) && time_idx <= length(texs[layer_idx])
            ny, nx = size(texs[layer_idx][time_idx])
            break
        end
    end

    # If no simulation has data for this layer/time, return empty matrices
    if ny == 0 || nx == 0
        # Get dimensions from any available layer
        for texs in all_texs
            if !isempty(texs) && !isempty(texs[1])
                ny, nx = size(texs[1][1])
                break
            end
        end
        return zeros(Float64, ny, nx), zeros(Int, ny, nx)
    end

    # Count how many simulations have CO2 at each location
    n_filled = zeros(Int, ny, nx)
    for (sim_idx, texs) in enumerate(all_texs)
        # Check if this simulation has data for this layer and time
        if layer_idx <= length(texs) && time_idx <= length(texs[layer_idx])
            state = texs[layer_idx][time_idx]
            # CO2 is present if state > 1 (1 = empty)
            n_filled .+= (state .> 1)
        end
    end

    # Compute fraction (0.0 to 1.0)
    fraction_filled = Float64.(n_filled) ./ n_simulations

    return fraction_filled, n_filled
end

"""
    compute_ensemble_statistics(all_texs)

Compute fraction of simulations with CO2 present across all layers and time points.

# Arguments
- `all_texs::Vector{Vector{Vector}}`: Textures from `run_ensemble_simulations`

# Returns
- `fraction_filled::Vector{Vector{Matrix{Float64}}}`: Fraction (0.0-1.0) for each layer and time
- `n_filled::Vector{Vector{Matrix{Int}}}`: Count of simulations with CO2 for each layer and time
"""
function compute_ensemble_statistics(all_texs::Vector)
    # Find maximum number of layers and times across all simulations
    n_layers = maximum(length.(all_texs))
    n_times = maximum(sim -> maximum(length.(sim)), all_texs)

    fraction_filled_all = []
    n_filled_all = []

    for layer_idx in 1:n_layers
        layer_fractions = []
        layer_counts = []

        for time_idx in 1:n_times
            fraction_filled, n_filled = compute_ensemble_statistics(all_texs, layer_idx, time_idx)
            push!(layer_fractions, fraction_filled)
            push!(layer_counts, n_filled)
        end

        push!(fraction_filled_all, layer_fractions)
        push!(n_filled_all, layer_counts)
    end

    return fraction_filled_all, n_filled_all
end
