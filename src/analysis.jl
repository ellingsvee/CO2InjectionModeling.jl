using SurfaceWaterIntegratedModeling: SpillEvent, WeatherEvent, inflow_at, numtraps,
    trap_states_at_timepoints

export LayerSnapshot, generate_layer_snapshot, generate_layer_snapshots
export total_to_next_layer, total_passthrough, total_upward_leakage
export compute_total_injected, compute_total_drained, compute_total_passthrough
export MultiLayerSnapshot, generate_multi_layer_snapshot, generate_multi_layer_snapshots
export print_summary

"""
A snapshot of a layer's CO2 state at a specific point in time, computed from the
fill sequence and leakage state. Designed to support mass balance analysis.
"""
struct LayerSnapshot
    layer_idx::Int
    layer_name::String
    timestamp::Float64

    # Per-trap state (indexed by trap ID, in SWIM volume units)
    trap_volumes::Vector{Float64}   # Current stored volume, corrected for residual drainage
    trap_filled::Vector{Bool}       # Structurally full (at spillpoint)
    trap_leaking::Vector{Bool}      # Has reached leakage threshold (edge = 0)
    trap_draining::Vector{Bool}     # Currently experiencing residual drainage
    trap_leakage_height::Vector{Float64}  # Per-trap CO2 column height threshold for leakage (Inf = cannot leak)

    # Mass balance quantities (all in SWIM volume units)
    total_injected::Float64         # Cumulative CO2 injected into this layer up to timestamp
    total_stored::Float64           # CO2 currently held in traps (= sum of trap_volumes)
    total_drained::Float64          # CO2 that drained out via residual leakage
    total_passthrough::Float64      # CO2 that flowed through leaking traps (fast path)
end

"""
Total CO2 that has left this layer up to `s.timestamp` via any path
(lateral domain spillage + upward leakage). Equal to `total_injected - total_stored`.
"""
total_to_next_layer(s::LayerSnapshot) = s.total_injected - s.total_stored

"""
    total_passthrough(s::LayerSnapshot) -> Float64

CO2 that flowed directly through leaking traps to the next layer (fast path),
as opposed to the slow residual drainage path.
"""
total_passthrough(s::LayerSnapshot) = s.total_passthrough

"""
    total_upward_leakage(s::LayerSnapshot) -> Float64

Total CO2 that leaked upward to the next layer via caprock failure.
Equals residual drainage + passthrough. Does NOT include lateral domain spillage.
"""
total_upward_leakage(s::LayerSnapshot) = s.total_drained + s.total_passthrough


"""
Integrate the injection rate over all weather events up to time `t`.
Returns total CO2 injected in SWIM volume units.
"""
function compute_total_injected(
    weather_events::Vector{WeatherEvent},
    t::Float64,
)::Float64
    total = 0.0
    for i in eachindex(weather_events)
        t_start = weather_events[i].timestamp
        t_start >= t && break

        t_end = (i == lastindex(weather_events)) ? t : min(weather_events[i+1].timestamp, t)
        duration = t_end - t_start
        duration <= 0 && continue

        rr = weather_events[i].rain_rate

        @assert rr isa Matrix "rain_rate must be a matrix of per-cell rates"
        total += sum(rr) * duration
    end
    return total
end

"""
Compute the total CO2 that flowed through leaking traps (passthrough) up to time `t`.
This is the "fast path": inflow to a leaking trap exits immediately upward to the next layer.
"""
function compute_total_passthrough(
    t::Float64,
    seq::Vector{SpillEvent},
    leakage_state::LeakageState,
    n_traps::Int
)::Float64
    total = 0.0
    for i in eachindex(seq)
        t_start = seq[i].timestamp
        t_start >= t && break
        t_end = (i == lastindex(seq)) ? t : min(seq[i+1].timestamp, t)
        duration = t_end - t_start
        duration <= 0 && continue
        inflows = inflow_at(seq, i)
        for trap in 1:n_traps
            !leakage_state.leaking[trap] && continue
            leakage_state.leakage_start_time[trap] > t_start && continue
            total += max(0.0, inflows[trap]) * duration
        end
    end
    return total
end

"""
Compute the total CO2 that has drained from traps via residual leakage up to time `t`.
This is the "slow path": CO2 that was stored in a trap and gradually drains after
the leakage threshold is crossed.
"""
function compute_total_drained(t::Float64, leakage_state::LeakageState)::Float64
    total = 0.0
    for trap in eachindex(leakage_state.draining)
        !leakage_state.draining[trap] && continue
        t < leakage_state.leakage_start_time[trap] && continue

        current_vol = compute_volume_with_drainage(trap, t, leakage_state)
        isnothing(current_vol) && continue

        drained = leakage_state.initial_volume_at_leak[trap] - current_vol
        total += max(0.0, drained)
    end
    return total
end

"""
Get per-trap volumes at time `t`, correcting SWIM's output for residual drainage in
draining traps. Returns four vectors indexed by trap ID.
"""
function _compute_trap_volumes_at_time(
    t::Float64,
    seq::Vector{SpillEvent},
    leakage_state::LeakageState,
    tstruct,
    tstate
)
    n = numtraps(tstruct)
    # tstates = trap_states_at_timepoints(tstruct, seq, [t]; verbose=false)
    swim_filled, swim_volumes, _ = tstate

    volumes = copy(swim_volumes)
    leaking = copy(leakage_state.leaking)
    draining = copy(leakage_state.draining)

    # Override SWIM volumes for traps that have started draining by time t.
    # SWIM keeps the volume constant at the value set when leakage was detected;
    # the true volume decreases linearly over residual_leakage_time.
    for trap in 1:n
        !draining[trap] && continue
        t < leakage_state.leakage_start_time[trap] && continue

        corrected = compute_volume_with_drainage(trap, t, leakage_state)
        isnothing(corrected) && continue
        volumes[trap] = corrected

        # If the trap hasn't started draining yet from SWIM's perspective
        # (leakage_start_time is in future relative to SWIM's seq), clear the flag.
        # (This shouldn't happen in practice but guards against edge cases.)
    end

    # Leaking/draining state: mask out traps whose leakage hadn't started yet at t
    for trap in 1:n
        if leaking[trap] && t < leakage_state.leakage_start_time[trap]
            leaking[trap] = false
        end
        if draining[trap] && t < leakage_state.leakage_start_time[trap]
            draining[trap] = false
        end
    end

    return volumes, swim_filled, leaking, draining
end


"""
Compute a [`LayerSnapshot`](@ref) for `layer` at time `t`.

# Arguments
- `layer`          : the analyzed layer (contains `trap_structure`)
- `layer_idx`      : index of this layer in the multi-layer stack (1 = deepest)
- `seq`            : fill-sequence produced by `fill_sequence_with_leakage`
- `leakage_state`  : leakage state produced by `fill_sequence_with_leakage`
- `weather_events` : the `WeatherEvent` vector passed to `fill_sequence_with_leakage`
- `t`              : the time at which to compute the snapshot (must be ≥ seq[1].timestamp)
"""
function generate_layer_snapshot(
    layer::Layer,
    layer_idx::Int,
    seq::Vector{SpillEvent},
    leakage_state::LeakageState,
    weather_events::Vector{WeatherEvent},
    t::Float64,
    tstate
)::LayerSnapshot

    tstruct = layer.trap_structure

    volumes, filled, leaking, draining =
        _compute_trap_volumes_at_time(t, seq, leakage_state, tstruct, tstate)

    total_stored = sum(volumes)
    total_inj = compute_total_injected(weather_events, t)
    total_drained_ = compute_total_drained(t, leakage_state)
    total_passthrough_ = compute_total_passthrough(t, seq, leakage_state, numtraps(tstruct))

    return LayerSnapshot(
        layer_idx,
        layer.name,
        t,
        volumes,
        filled,
        leaking,
        draining,
        leakage_state.leakage_height,
        total_inj,
        total_stored,
        total_drained_,
        total_passthrough_
    )
end

"""
    generate_layer_snapshots(layer, layer_idx, seq, leakage_state,
                             weather_events, timepoints) -> Vector{LayerSnapshot}

Compute [`LayerSnapshot`](@ref)s for `layer` at every time in `timepoints`.

SWIM trap states are evaluated in bulk via `trap_states_at_timepoints` for
efficiency, making this much faster than calling [`generate_layer_snapshot`](@ref)
in a loop.  This is the preferred method when multiple snapshots are needed.

# Arguments
- `layer`: [`Layer`](@ref)
- `layer_idx`: Index of this layer in the multi-layer stack (1 = deepest)
- `seq`: Fill sequence from [`fill_sequence_with_leakage`](@ref)
- `leakage_state`: Leakage state from [`fill_sequence_with_leakage`](@ref)
- `weather_events`: Effective weather events for this layer
- `timepoints`: Sorted vector of query times
"""
function generate_layer_snapshots(
    layer::Layer,
    layer_idx::Int,
    seq::Vector{SpillEvent},
    leakage_state::LeakageState,
    weather_events::Vector{WeatherEvent},
    timepoints::Vector{Float64}
)::Vector{LayerSnapshot}

    tstruct = layer.trap_structure
    tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)
    return [generate_layer_snapshot(layer, layer_idx, seq, leakage_state, weather_events, t, tstates[i])
            for (i, t) in enumerate(timepoints)]
end


"""
A snapshot of the entire multi-layer CO2 stack at a specific point in time.

# Fields
- `timestamp`              : simulation time
- `layers`                 : per-layer snapshots (index 1 = deepest / injection layer)
- `total_injected`         : cumulative CO2 injected into the system (= `layers[1].total_injected`)
- `total_stored`           : CO2 currently held in all traps across all layers
- `total_surface_leakage`  : CO2 that has escaped from the top layer

# Key mass-conservation invariant (closed BC, injection only in layer 1):
```
total_injected ≈ total_stored + total_surface_leakage
```
"""
struct MultiLayerSnapshot
    timestamp::Float64
    layers::Vector{LayerSnapshot}   # index 1 = deepest (injection) layer
    total_injected::Float64         # = layers[1].total_injected
    total_stored::Float64           # sum(s.total_stored for s in layers)
    total_surface_leakage::Float64  # total_to_next_layer(layers[end])
end

"""
    generate_multi_layer_snapshot(layers, seqs, leakage_states,
                                  weather_events_per_layer, t) -> MultiLayerSnapshot

Compute a [`MultiLayerSnapshot`](@ref) for the entire layer stack at time `t`.

# Arguments
- `layers`: `Vector{Layer}` (deepest first)
- `seqs`: Fill sequences from [`fill_layers`](@ref)
- `leakage_states`: Leakage states from [`fill_layers`](@ref)
- `weather_events_per_layer`: Weather events per layer from [`fill_layers`](@ref)
- `t`: Query time

See also [`generate_multi_layer_snapshots`](@ref) for computing many snapshots
efficiently.
"""
function generate_multi_layer_snapshot(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    leakage_states::Vector{LeakageState},
    weather_events_per_layer::Vector{Vector{WeatherEvent}},
    t::Float64
)::MultiLayerSnapshot
    n = length(layers)
    layer_snaps = Vector{LayerSnapshot}(undef, n)
    for i in 1:n
        tstruct = layers[i].trap_structure
        tstate = trap_states_at_timepoints(tstruct, seqs[i], [t]; verbose=false)[1]
        layer_snaps[i] = generate_layer_snapshot(
            layers[i], i, seqs[i], leakage_states[i], weather_events_per_layer[i], t, tstate
        )
    end

    total_inj = layer_snaps[1].total_injected
    total_stored = sum(s.total_stored for s in layer_snaps)
    total_surface_leakage = total_to_next_layer(layer_snaps[end])

    return MultiLayerSnapshot(t, layer_snaps, total_inj, total_stored, total_surface_leakage)
end

"""
    generate_multi_layer_snapshots(layers, seqs, leakage_states,
                                   weather_events_per_layer, timepoints)
        -> Vector{MultiLayerSnapshot}

Compute [`MultiLayerSnapshot`](@ref)s for the entire layer stack at each time
in `timepoints`.

Trap states are pre-computed for all timepoints in bulk (per layer) for
efficiency.  This is the standard entry point for post-processing after
[`fill_layers`](@ref).

# Arguments
- `layers`: `Vector{Layer}` (deepest first)
- `seqs`, `leakage_states`, `weather_events_per_layer`: outputs of [`fill_layers`](@ref)
- `timepoints`: Sorted vector of query times

# Example
```julia
timepoints = collect(range(0.0, 15.0, length=30))
multi_snaps = generate_multi_layer_snapshots(
    layers, seqs, leakage_states, weather_events_per_layer, timepoints)
print_summary(multi_snaps[end])
```
"""
function generate_multi_layer_snapshots(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    leakage_states::Vector{LeakageState},
    weather_events_per_layer::Vector{Vector{WeatherEvent}},
    timepoints::Vector{Float64}
)::Vector{MultiLayerSnapshot}
    n = length(layers)
    # Pre-compute trap states for each layer over all timepoints
    layer_snap_vectors = Vector{Vector{LayerSnapshot}}(undef, n)
    for i in 1:n
        layer_snap_vectors[i] = generate_layer_snapshots(
            layers[i], i, seqs[i], leakage_states[i], weather_events_per_layer[i], timepoints
        )
    end

    return [
        begin
            snaps_at_t = [layer_snap_vectors[i][ti] for i in 1:n]
            total_inj = snaps_at_t[1].total_injected
            total_stored = sum(s.total_stored for s in snaps_at_t)
            total_surface_leakage = total_to_next_layer(snaps_at_t[end])
            MultiLayerSnapshot(t, snaps_at_t, total_inj, total_stored, total_surface_leakage)
        end
        for (ti, t) in enumerate(timepoints)
    ]
end


"""
    print_summary([io,] snap::MultiLayerSnapshot)

Print a compact human-readable mass-balance table for `snap` to `io`
(default: `stdout`).

Columns: layer name, cumulative injected, currently stored, residually drained,
passthrough, and volume leaving to the next layer.  A mass-imbalance line is
printed at the bottom.

# Example
```julia
print_summary(multi_snaps[end])
```
"""
function print_summary(io::IO, snap::MultiLayerSnapshot)
    @printf(io, "MultiLayerSnapshot  t = %.4g\n", snap.timestamp)
    @printf(io, "  %-6s  %12s  %12s  %12s  %12s  %12s\n",
        "Layer", "Injected", "Stored", "Drained", "Passthru", "→NextLayer")
    @printf(io, "  %s\n", "-"^74)
    for s in snap.layers
        @printf(io, "  %-6s  %12.4g  %12.4g  %12.4g  %12.4g  %12.4g\n",
            s.layer_name,
            s.total_injected,
            s.total_stored,
            s.total_drained,
            s.total_passthrough,
            total_to_next_layer(s))
    end
    @printf(io, "  %s\n", "-"^74)
    imbalance = snap.total_injected - snap.total_stored - snap.total_surface_leakage
    @printf(io, "  %-6s  %12.4g  %12.4g  %12s  %12s  %12.4g\n",
        "TOTAL", snap.total_injected, snap.total_stored, "", "", snap.total_surface_leakage)
    @printf(io, "  Mass imbalance: %.4g\n", imbalance)
end

print_summary(snap::MultiLayerSnapshot) = print_summary(stdout, snap)


"""
    print_summary([io,] snap, reservoir_properties, domain)

Print a mass-balance table in physical mass units (Mt) using per-layer CO2
density and pore-volume scaling.  This is the correct way to display mass
balance when layers have different CO2 densities.

The mass imbalance should be zero (within numerical tolerance) regardless of
density differences between layers.
"""
function print_summary(
    io::IO,
    snap::MultiLayerSnapshot,
    reservoir_properties::Vector{ReservoirProperties},
    domain::Domain3D,
)
    @printf(io, "MultiLayerSnapshot  t = %.4g  (mass units: Mt)\n", snap.timestamp)
    @printf(io, "  %-6s  %12s  %12s  %12s  %12s  %12s\n",
        "Layer", "Injected", "Stored", "Drained", "Passthru", "→NextLayer")
    @printf(io, "  %s\n", "-"^74)

    total_injected_mt = 0.0
    total_stored_mt = 0.0
    total_surface_leakage_mt = 0.0

    for (i, s) in enumerate(snap.layers)
        rp = reservoir_properties[i]
        scale = full_volume_to_rock_volume_scaling(rp) * unit_volume_to_physical_scaling(domain)
        to_mt = scale * rp.co2_density / 1e9

        inj_mt = s.total_injected * to_mt
        stored_mt = s.total_stored * to_mt
        drained_mt = s.total_drained * to_mt
        passthru_mt = s.total_passthrough * to_mt
        next_mt = total_to_next_layer(s) * to_mt

        @printf(io, "  %-6s  %12.4g  %12.4g  %12.4g  %12.4g  %12.4g\n",
            s.layer_name, inj_mt, stored_mt, drained_mt, passthru_mt, next_mt)

        if i == 1
            total_injected_mt = inj_mt
        end
        total_stored_mt += stored_mt
        if i == length(snap.layers)
            total_surface_leakage_mt = next_mt
        end
    end

    @printf(io, "  %s\n", "-"^74)
    imbalance = total_injected_mt - total_stored_mt - total_surface_leakage_mt
    @printf(io, "  %-6s  %12.4g  %12.4g  %12s  %12s  %12.4g\n",
        "TOTAL", total_injected_mt, total_stored_mt, "", "", total_surface_leakage_mt)
    @printf(io, "  Mass imbalance: %.4g Mt\n", imbalance)
end

function print_summary(
    snap::MultiLayerSnapshot,
    reservoir_properties::Vector{ReservoirProperties},
    domain::Domain3D,
)
    print_summary(stdout, snap, reservoir_properties, domain)
end
