using SurfaceWaterIntegratedModeling

export SimulationSummary, generate_simulation_summary

"""
Structure to hold summary statistics from a simulation run.

# Fields
- `timepoints::Vector{Float64}`: Vector of timestamps from the snapshots
- `total_co2_volumes::Vector{Float64}`: Total CO2 volume at each timepoint
- `layer_co2_volumes::Matrix{Float64}`: CO2 volume in each layer at each timepoint (timepoints × layers)
- `trap_co2_volumes::Vector{Matrix{Float64}}`: CO2 volume in each trap for each layer at each timepoint
  (vector of length = num_layers, each matrix is timepoints × num_traps_in_layer)
- `trap_co2_percentages::Vector{Matrix{Float64}}`: CO2 volume in each trap as percentage of total volume
  (same structure as trap_co2_volumes)
- `num_layers::Int`: Number of layers in the simulation
- `num_traps_per_layer::Vector{Int}`: Number of traps in each layer
"""
struct SimulationSummary
    timepoints::Vector{Float64}
    total_co2_volumes::Vector{Float64}
    layer_co2_volumes::Matrix{Float64}
    trap_co2_volumes::Vector{Matrix{Float64}}
    trap_co2_percentages::Vector{Matrix{Float64}}
    num_layers::Int
    num_traps_per_layer::Vector{Int}
end

"""
    generate_simulation_summary(snapshots::Vector{SimulationSnapshot}, layers::Vector{Layer}, seqs::Vector)

Generate a comprehensive summary from simulation snapshots.

# Arguments
- `snapshots::Vector{SimulationSnapshot}`: Vector of snapshots from fill_layers
- `layers::Vector{Layer}`: Vector of layer structures used in the simulation
- `seqs::Vector`: Vector of SpillEvent sequences from fill_layers (one sequence per layer)

# Returns
- `SimulationSummary`: Structure containing time series of volumes and percentages

# Example
```julia
seqs, leakage_states, snapshots = fill_layers(
    layers, domain, reservoir_props, injection_events;
    num_snapshots=14, start_time=0.0, end_time=15.0
)

summary = generate_simulation_summary(snapshots, layers, seqs)

# Access the data
println("Timepoints: ", summary.timepoints)
println("Total CO2 at each time: ", summary.total_co2_volumes)
println("Layer 1 volumes: ", summary.layer_co2_volumes[:, 1])
println("Layer 1, Trap 3 volumes: ", summary.trap_co2_volumes[1][:, 3])
println("Layer 1, Trap 3 percentages: ", summary.trap_co2_percentages[1][:, 3])
```
"""
function generate_simulation_summary(
    snapshots::Vector{SimulationSnapshot},
    layers::Vector{Layer},
    seqs::Vector
)::SimulationSummary

    num_snapshots = length(snapshots)
    num_layers = length(layers)

    # Get number of traps in each layer
    num_traps_per_layer = [numtraps(layer.trap_structure) for layer in layers]

    # Initialize arrays
    timepoints = zeros(Float64, num_snapshots)
    total_co2_volumes = zeros(Float64, num_snapshots)
    layer_co2_volumes = zeros(Float64, num_snapshots, num_layers)

    # Initialize trap volumes - one matrix per layer
    trap_co2_volumes = [zeros(Float64, num_snapshots, num_traps_per_layer[i])
                        for i in 1:num_layers]
    trap_co2_percentages = [zeros(Float64, num_snapshots, num_traps_per_layer[i])
                           for i in 1:num_layers]

    # Process each snapshot
    for (snap_idx, snapshot) in enumerate(snapshots)
        timepoints[snap_idx] = snapshot.timestamp
        total_co2_volumes[snap_idx] = snapshot.total_co2_volume

        # Process each layer in this snapshot
        for (layer_idx, layer_snapshot) in enumerate(snapshot.layer_snapshots)
            layer_co2_volumes[snap_idx, layer_idx] = layer_snapshot.co2_volume

            # Find the event index in the sequence that corresponds to this snapshot timestamp
            # We need to find which event in seqs[layer_idx] matches this timestamp
            seq = seqs[layer_idx]
            event_idx = findfirst(e -> e.timestamp == snapshot.timestamp, seq)

            # Extract volumes for each trap using amount_at
            # This handles both Vector{FilledAmount} and Vector{IncrementalUpdate{FilledAmount}}
            if !isnothing(event_idx)
                amounts = amount_at(seq, event_idx)
            else
                # If exact match not found, use the last event before or at this timestamp
                event_idx = findlast(e -> e.timestamp <= snapshot.timestamp, seq)
                if !isnothing(event_idx)
                    amounts = amount_at(seq, event_idx)
                else
                    # No events yet, use zeros
                    amounts = [SurfaceWaterIntegratedModeling.FilledAmount(0.0, 0.0)
                              for _ in 1:num_traps_per_layer[layer_idx]]
                end
            end

            for trap_idx in 1:num_traps_per_layer[layer_idx]
                trap_volume = amounts[trap_idx].amount
                trap_co2_volumes[layer_idx][snap_idx, trap_idx] = trap_volume

                # Calculate percentage of total volume
                if snapshot.total_co2_volume > 0.0
                    trap_co2_percentages[layer_idx][snap_idx, trap_idx] =
                        (trap_volume / snapshot.total_co2_volume) * 100.0
                else
                    trap_co2_percentages[layer_idx][snap_idx, trap_idx] = 0.0
                end
            end
        end
    end

    return SimulationSummary(
        timepoints,
        total_co2_volumes,
        layer_co2_volumes,
        trap_co2_volumes,
        trap_co2_percentages,
        num_layers,
        num_traps_per_layer
    )
end
