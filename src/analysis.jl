# Some structs for "snapshots" of the simulation state at a given time, used for analysis and visualization

# Note the SurfaceWaterIntegratedModeling.trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false) is very usefull for getting trap states at spesific time-points.

struct LayerSnapshot
    # Basic info
    layer_idx::Int
    layer_name::String
    timestamp::Float64

    # Trap states
    num_traps::Int
    num_filled_traps::Int
    # Maybe something related to the leakage and drainage states?

    # Volumes
    total_injected_volume::Float64 # In SWIM units
    total_stored_volume::Float64 # In SWIM units
    total_drained_volume::Float64 # In SWIM units

    # Maybe something more as well
end


function generate_layer_snapshots(
    layer::Layer,
    seqs::Vector{SpillEvent},
    leakage_state::LeakageState,
    timepoints::Vector{Float64}
)
    layer_snapshots = Vector{LayerSnapshot}(undef, length(timepoints))
    trap_states = SurfaceWaterIntegratedModeling.trap_states_at_timepoints(layer.trap_structure, seqs, timepoints; verbose=false)[1]


    for (i, t) in enumerate(timepoints)
        # Get trap states at this timepoint
        filled, volumes, _ = trap_states[i]

        # TODO: Not finished
    end

    return layer_snapshots
end