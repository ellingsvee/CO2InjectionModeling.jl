module CO2RInterface

using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using Interpolations

export SimulatorState, setup_simulator, run_simulation

"""
Struct holding all data required to run a simulation.
This keeps the state clean and easy to pass from R.
"""
struct SimulatorState
    layers
    domain
    reservoir_props
    base_topography   # trap_structure.topography from the bottom layer
end


"""
setup_simulator()

Creates all geological structures, lithology, domain, layers,
and reservoir properties. This is the heavy initialization.

Returns a SimulatorState object.
"""
function setup_simulator(; topography_scale = 1.0)

    # Load topography and geological structures
    topography = load_sleipner_topography()
    domain = create_domain_from_topography(topography, topography_scale)
    lithology = reconstruct_3d_lithology(topography, domain)
    layers = analyze_base_surfaces(topography)

    reservoir_props = ReservoirProperties()

    # Most users inject into the bottom layer → layer 1
    base_topography = layers[1].trap_structure.topography

    return SimulatorState(layers, domain, reservoir_props, base_topography)
end



"""
run_simulation(sim, injection_cells, injection_values; num_snapshots, start_time, end_time)

Inputs:
- sim::SimulatorState  – Output from setup_simulator()
- injection_cells      – Vector of (i,j) coordinates or a 2-column Matrix
- injection_values     – Vector of numeric injection rates (same length as injection_cells)

Keyword args:
- num_snapshots = 10
- start_time = 0.0
- end_time = 10.0

Returns:
- seqs, leakage_states, snapshots

This is the function intended to be called from R.
"""
function run_simulation(
        sim::SimulatorState,
        injection_cells,
        injection_values;
        num_snapshots = 10,
        start_time = 0.0,
        end_time = 10.0,
        verbose = false
    )

    layers = sim.layers
    domain = sim.domain
    reservoir_props = sim.reservoir_props

    # Prepare injection arrays for each layer (only layer 1 gets injections)
    base_topo = sim.base_topography
    injection_rate_bottom = fill(0.0, size(base_topo))

    # Accept either:
    #  - injection_cells = [(i,j), ...] or
    #  - injection_cells = Matrix{Int}(N×2)
    if isa(injection_cells, AbstractMatrix)
        @assert size(injection_cells, 2) == 2 "Matrix must be N×2 (i,j)"
        for k in 1:size(injection_cells,1)
            i, j = injection_cells[k, :]
            injection_rate_bottom[i, j] = injection_values[k]
        end
    else
        @assert length(injection_cells) == length(injection_values)
        for (idx,(i,j)) in enumerate(injection_cells)
            injection_rate_bottom[i, j] = injection_values[idx]
        end
    end

    # Construct injection events
    injection_event_bottom_layer = [
        InjectionEvent(start_time, injection_rate_bottom)
    ]

    # All other layers → zero injection
    zero_events = [
        InjectionEvent(start_time, zeros(size(base_topo)))
    ]

    injection_events = Vector{Vector{InjectionEvent}}(undef, length(layers))
    injection_events[1] = injection_event_bottom_layer
    for i in 2:length(layers)
        injection_events[i] = zero_events
    end

    # Run the model
    seqs, leakage_states, snapshots = fill_layers(
        layers,
        domain,
        reservoir_props,
        injection_events;
        num_snapshots = num_snapshots,
        start_time = start_time,
        end_time = end_time,
        verbose = verbose
    )

    return seqs, leakage_states, snapshots
end

end # module
