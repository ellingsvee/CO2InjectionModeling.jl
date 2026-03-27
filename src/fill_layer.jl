using SurfaceWaterIntegratedModeling
using Statistics: mean, std
export fill_sequence_with_leakage, InjectionEvent, get_effective_leakage_cap


"""
    fill_sequence_with_leakage(tstruct, reservoir_properties, weather_events;
                               time_slack=0.0, verbose=false)
        -> (seq::Vector{SpillEvent}, leakage_state::LeakageState)

Simulate CO2 filling for a single layer represented by `tstruct`, driven by
`weather_events` (injection rates in SWIM volume units).  Leakage through the
caprock is detected and tracked: when the CO2 column height in a trap exceeds
the threshold in `reservoir_properties`, that trap is flagged as leaking and
its CO2 begins to drain residually.

Normally called indirectly through [`fill_layers`](@ref).

# Arguments
- `tstruct`: SWIM `TrapStructure` for this layer (from [`analyze_base_surfaces`](@ref))
- `reservoir_properties`: [`ReservoirProperties`](@ref) for this layer
- `weather_events`: `Vector{WeatherEvent}` giving injection rates over time
  (produced by [`convert_injection_event_to_weather_event`](@ref))
- `time_slack`: Small time offset added to SWIM events (default `0.0`)
- `verbose`: Print progress information (default `false`)

# Returns
- `seq`: `Vector{SpillEvent}` describing the trap-filling sequence (SWIM output)
- `leakage_state`: [`LeakageState`](@ref) recording which traps leaked, when,
  and how much
"""
function fill_sequence_with_leakage(tstruct::TrapStructure{<:Real},
    reservoir_properties::ReservoirProperties,
    weather_events::Vector{WeatherEvent};
    time_slack::Float64=0.0,
    verbose::Bool=false)


    @assert !isempty(weather_events)

    num_traps = numtraps(tstruct)
    if num_traps == 0
        # No traps - return empty results
        empty_leakage = LeakageState(
            Bool[],  # leaking
            Bool[],  # draining
            Float64[], Float64[], LeakageRecord[],
            Float64[],  # leakage_height (per-trap vector, empty for no traps)
            Float64[],  # initial_volume_at_leak
            reservoir_properties.sand_residual_co2_saturation /
                (1.0 - reservoir_properties.sand_irreducible_water_saturation),
            reservoir_properties.residual_leakage_time,
            Float64[],  # cumulative_no_inflow_time
            Float64[],  # volume_at_last_state_change
            Float64[],  # time_of_last_state_change
            Bool[]      # has_inflow
        )
        return Vector{SpillEvent}(), empty_leakage
    end

    # Initialize infiltration map from user input
    infiltration = zeros(size(tstruct.topography))

    # Compute tables to support computation of trap water volume as function of water level
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # Set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), num_traps)
    sgraph = SurfaceWaterIntegratedModeling.compute_complete_spillgraph(tstruct, filled_traps)

    # Start with empty sequence
    seq = Vector{SpillEvent}()

    # Initialize leakage state
    leakage_state = initialize_leakage_state(
        tstruct, z_vol_tables, reservoir_properties,
        reservoir_properties.sand_residual_co2_saturation,
        reservoir_properties.residual_leakage_time
    )

    if verbose
        finite_heights = leakage_state.leakage_height[isfinite.(leakage_state.leakage_height)]
        if isempty(finite_heights)
            println("Leakage height: sealed (all Inf)")
        else
            mean_height = mean(finite_heights)
            std_height = std(finite_heights)
            println("Leakage height: mean=$(round(mean_height, digits=2)) m, std=$(round(std_height, digits=2)) m")
        end
    end

    # Compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp

        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # Compute inflow/runoff/infiltration rates corresponding to the fill graph and new rain rate
        rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, we.rain_rate, infiltration, tstruct, verbose)

        # Update dynamic equilibrium state for leaking traps at weather event boundary.
        # This detects when injection rate changes (e.g., goes to 0), updating has_inflow.
        for trap in 1:num_traps
            leakage_state.leaking[trap] || continue
            new_inflow = SurfaceWaterIntegratedModeling.getinflow(rateinfo, trap) > 0
            update_leaking_trap_inflow_state!(leakage_state, trap, new_inflow, cur_time)
        end

        # Compute initial time estimates for when a trap become filled, or split into subtraps
        changetimeest = SurfaceWaterIntegratedModeling._set_initial_changetime_estimates(rateinfo, cur_amounts,
            cur_time, filled_traps,
            tstruct)

        # Compute initial leakage time estimates
        leakage_time_est = get_initial_leakage_time_estimates(
            cur_amounts, cur_time, rateinfo, leakage_state, filled_traps, tstruct
        )

        # Register the start of this weather event as a new, fully computed, spill event
        push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
            copy(rateinfo.trap_inflow), copy(we.rain_rate),
            copy(rateinfo.runoff)))

        # Will add new events to `seq`. `sgraph`, `rateinfo`, `changetimeest`,
        # `filled_traps`, `cur_amounts`, and `leakage_state` are modified in the process
        _fill_sequence_for_weather_event_with_leakage!(
            seq, sgraph, rateinfo, changetimeest, leakage_time_est,
            filled_traps, cur_amounts, z_vol_tables,
            tstruct, infiltration, end_time, time_slack,
            leakage_state, verbose
        )
    end

    return seq, leakage_state
end


"""
    _fill_sequence_for_weather_event_with_leakage!(...)

Modified fill sequence that handles both fill/empty events and leakage events.
"""
function _fill_sequence_for_weather_event_with_leakage!(
    seq, sgraph, rateinfo, changetimeest, leakage_time_est,
    filled_traps, cur_amounts, z_vol_tables,
    tstruct, infiltration, endtime, time_slack,
    leakage_state::LeakageState, verbose
)
    cur_time = cur_amounts[1].time
    num_traps = numtraps(tstruct)

    fill_updates = Vector{IncrementalUpdate{Bool}}()
    graph_updates = Vector{IncrementalUpdate{Int}}()

    count = 0
    while cur_time < endtime
        verbose && (mod(count += 1, 10) == 0) && println("Fill sequence iteration: ", count)

        # Find next fill/empty event (from SWIM)
        next_fill_time, fill_updates = SurfaceWaterIntegratedModeling._identify_next_status_change!(
            changetimeest, cur_amounts, rateinfo,
            filled_traps, tstruct, z_vol_tables,
            cur_time, endtime
        )

        # Find next leakage event
        next_leak_time, leak_trap = find_next_leakage_event(leakage_time_est)

        # Determine which event comes first
        if next_leak_time < next_fill_time && next_leak_time <= endtime && leak_trap > 0
            # LEAKAGE EVENT occurs first
            cur_time = next_leak_time

            verbose && println("Leakage event at time $(cur_time) in trap $(leak_trap)")

            # Mark trap as leaking and record start time.
            leakage_state.leaking[leak_trap] = true
            leakage_state.leakage_start_time[leak_trap] = cur_time

            # Initialize dynamic equilibrium: check if trap has inflow at leakage onset.
            # The trap just reached threshold, so it has inflow (otherwise it wouldn't fill).
            inflow_at_leak = SurfaceWaterIntegratedModeling.getinflow(rateinfo, leak_trap)
            leakage_state.has_inflow[leak_trap] = inflow_at_leak > 0
            leakage_state.time_of_last_state_change[leak_trap] = cur_time
            # volume_at_last_state_change will be set below after leakage_vol is computed

            # Record the leakage for upstream layer
            leakage_location = find_leakage_location(leak_trap, tstruct)
            push!(leakage_state.leakage_records, LeakageRecord(
                cur_time,
                leak_trap,
                leakage_location
            ))

            # Compute the actual volume at leakage time (before capping).
            leakage_vol = leakage_state.leakage_volume[leak_trap]
            actual_vol_at_leak = SurfaceWaterIntegratedModeling._compute_exact_fill(
                rateinfo, cur_amounts, leak_trap,
                filled_traps, tstruct, cur_time,
                z_vol_tables, false
            )
            initial_vol = max(leakage_vol, actual_vol_at_leak)
            leakage_state.initial_volume_at_leak[leak_trap] = initial_vol
            leakage_state.volume_at_last_state_change[leak_trap] = leakage_vol
            verbose && println("  Trap $leak_trap: leakage_vol=$(round(leakage_vol, digits=4)), initial_vol=$(round(initial_vol, digits=4))")

            # Mark as draining if there is actual volume to drain
            if initial_vol > 0
                leakage_state.draining[leak_trap] = true
            end

            # Mark all filled descendants as draining too. 
            descendants = get_all_descendants(tstruct, leak_trap)
            for desc_id in descendants
                if filled_traps[desc_id] && !leakage_state.draining[desc_id] && !leakage_state.leaking[desc_id]
                    leakage_state.draining[desc_id] = true
                    leakage_state.leakage_start_time[desc_id] = cur_time
                    # Get descendant's volume at leak time
                    desc_vol = cur_amounts[desc_id].amount
                    leakage_state.initial_volume_at_leak[desc_id] = desc_vol
                    verbose && println("  Descendant trap $(desc_id) marked as draining (vol=$(round(desc_vol, digits=2)))")
                end
            end

            # Cap the trap amount at leakage volume
            cur_amounts[leak_trap] = FilledAmount(leakage_vol, cur_time)

            # Mark only the leaking trap as filled with edge=0.
            # Ancestors keep their normal edges as they spill to this trap, which then leaks.
            filled_traps[leak_trap] = true

            verbose && println("  Trap $(leak_trap) marked as leaking, edge set to 0 (leak out of domain)")

            # Create fill update for the leaking trap only
            leak_fill_updates = [IncrementalUpdate{Bool}(leak_trap, true)]

            # Update spillgraph for the leaking trap
            graph_updates = SurfaceWaterIntegratedModeling.update_spillgraph!(sgraph, leak_fill_updates, tstruct)

            # Fix edges for all previously leaking traps (sibling cycle check may overwrite them)
            _fix_leaking_trap_edges!(sgraph, graph_updates, leakage_state)

            # Set only the leaking trap's edge to 0 (out of domain = leakage)
            sgraph.edges[leak_trap] = 0

            # Update flow information with the modified spillgraph
            setsavepoint!(rateinfo)
            SurfaceWaterIntegratedModeling._update_flow!(rateinfo, graph_updates, tstruct, sgraph)

            # Update amounts for all traps whose inflow rate changed due to the
            # spillgraph modification.  
            amount_updates = SurfaceWaterIntegratedModeling._update_affected_amounts(
                rateinfo, cur_amounts, filled_traps,
                tstruct, z_vol_tables, cur_time
            )

            # Amounts for unfilled descendants of the leaking trap.
            affected_set = Set(u.index for u in amount_updates)
            for desc_id in descendants
                if !filled_traps[desc_id] && !(desc_id in affected_set)
                    desc_vol = SurfaceWaterIntegratedModeling._compute_exact_fill(
                        rateinfo, cur_amounts, desc_id,
                        filled_traps, tstruct, cur_time,
                        z_vol_tables, true
                    )
                    if desc_vol > 0
                        push!(amount_updates, IncrementalUpdate(desc_id, FilledAmount(desc_vol, cur_time)))
                        verbose && println("  Materialized unfilled descendant $desc_id: vol=$(round(desc_vol, digits=4))")
                    end
                end
            end

            # Also cap amounts for any leaking traps that were updated by
            # _update_affected_amounts (e.g., pass-through parent traps).
            for i in 1:length(amount_updates)
                tix = amount_updates[i].index
                if leakage_state.leaking[tix]
                    cap_vol = get_effective_leakage_cap(leakage_state, tix)
                    if amount_updates[i].value.amount > cap_vol
                        amount_updates[i] = IncrementalUpdate(tix, FilledAmount(cap_vol, cur_time))
                    end
                end
            end

            # Add the amount update for the leaking trap itself
            push!(amount_updates, IncrementalUpdate(leak_trap, FilledAmount(leakage_vol, cur_time)))

            # Integrate amount changes into cur_amounts (must happen before
            # leakage/changetime estimate updates which read cur_amounts)
            SurfaceWaterIntegratedModeling._apply_updates!(cur_amounts, amount_updates)

            # Update dynamic equilibrium state for all leaking traps whose inflow changed
            for trap in 1:num_traps
                leakage_state.leaking[trap] || continue
                new_inflow = SurfaceWaterIntegratedModeling.getinflow(rateinfo, trap) > 0
                update_leaking_trap_inflow_state!(leakage_state, trap, new_inflow, cur_time)
            end

            # Update leakage time estimate to Inf for the leaking trap (already leaking)
            leakage_time_est[leak_trap] = LeakageTimeEstimate(leak_trap, Inf, Inf)

            # Update leakage time estimates for all traps whose inflow changed
            affected_traps = unique([u.index for u in getinflowupdates(rateinfo)])
            update_leakage_time_estimates!(
                leakage_time_est, affected_traps, cur_amounts, cur_time,
                rateinfo, leakage_state, filled_traps, tstruct
            )

            # Recompute SWIM's changetime estimates for all traps.
            for trap in 1:num_traps
                changetimeest[trap] = SurfaceWaterIntegratedModeling._compute_changetime_estimate(
                    trap, cur_amounts, cur_time, rateinfo, filled_traps, tstruct
                )
            end

            # Record this as a spill event
            push!(seq, SpillEvent(cur_time, amount_updates, leak_fill_updates,
                getinflowupdates(rateinfo), nothing,
                getrunoffupdates(rateinfo)))

        elseif next_fill_time <= endtime && !isempty(fill_updates)
            # FILL/EMPTY EVENT occurs first (standard SWIM logic)
            cur_time = next_fill_time

            for u in fill_updates
                filled_traps[u.index] = u.value
            end

            # Given changes in fill state, update spill graph
            graph_updates = SurfaceWaterIntegratedModeling.update_spillgraph!(sgraph, fill_updates, tstruct)

            # Fix edges for all leaking traps (sibling cycle check may overwrite them)
            _fix_leaking_trap_edges!(sgraph, graph_updates, leakage_state)

            # Given the updates to the spill graph, update flow information in `rateinfo`
            setsavepoint!(rateinfo)
            SurfaceWaterIntegratedModeling._update_flow!(rateinfo, graph_updates, tstruct, sgraph)

            # Update water amount in traps whose inflow rate is about to change, or that just filled
            amount_updates = SurfaceWaterIntegratedModeling._update_affected_amounts(
                rateinfo, cur_amounts, filled_traps,
                tstruct, z_vol_tables, cur_time
            )

            # For traps that just filled, set their amount to full capacity
            for tix in [u.index for u in fill_updates]
                if leakage_state.leaking[tix]
                    # Leaking trap - cap at leakage volume (or 0 if pass-through)
                    cap_vol = get_effective_leakage_cap(leakage_state, tix)
                    push!(amount_updates, IncrementalUpdate(tix, FilledAmount(cap_vol, cur_time)))
                else
                    # Normal trap - fill to capacity
                    fill_vol = tstruct.trapvolumes[tix] - tstruct.subvolumes[tix]
                    push!(amount_updates, IncrementalUpdate(tix, FilledAmount(fill_vol, cur_time)))

                    # Check if this newly filled trap feeds into a leaking or draining chain.
                    # A leaking trap includes pass-through traps (leakage_vol=0).
                    if !leakage_state.draining[tix]
                        spill_target = sgraph.edges[tix]
                        # Only check if spill_target is a valid trap index (not runoff/boundary)
                        if spill_target > 0 && spill_target <= num_traps &&
                           (leakage_state.leaking[spill_target] || leakage_state.draining[spill_target])
                            # This trap spills into a leaking/draining trap - it should also drain
                            leakage_state.draining[tix] = true
                            leakage_state.leakage_start_time[tix] = cur_time
                            leakage_state.initial_volume_at_leak[tix] = fill_vol
                            verbose && println("  Newly filled trap $(tix) marked as draining (feeds into leaking/draining trap $(spill_target))")
                        end
                    end
                end
            end

            # Also cap amounts for leaking traps that were updated by
            # _update_affected_amounts. This handles parent traps that are in "pass-through" mode.
            for i in 1:length(amount_updates)
                tix = amount_updates[i].index
                if leakage_state.leaking[tix]
                    cap_vol = get_effective_leakage_cap(leakage_state, tix)
                    if amount_updates[i].value.amount > cap_vol
                        amount_updates[i] = IncrementalUpdate(tix, FilledAmount(cap_vol, cur_time))
                    end
                end
            end

            # Integrate the changes into the continuously updated `cur_amounts` vector
            SurfaceWaterIntegratedModeling._apply_updates!(cur_amounts, amount_updates)

            # Update dynamic equilibrium state for all leaking traps whose inflow changed
            for trap in 1:num_traps
                leakage_state.leaking[trap] || continue
                new_inflow = SurfaceWaterIntegratedModeling.getinflow(rateinfo, trap) > 0
                update_leaking_trap_inflow_state!(leakage_state, trap, new_inflow, cur_time)
            end

            # Update leakage time estimates for affected traps
            affected_traps = unique([u.index for u in getinflowupdates(rateinfo)])
            update_leakage_time_estimates!(
                leakage_time_est, affected_traps, cur_amounts, cur_time,
                rateinfo, leakage_state, filled_traps, tstruct
            )

            # Add current state to result
            push!(seq, SpillEvent(cur_time, amount_updates, fill_updates,
                getinflowupdates(rateinfo), nothing,
                getrunoffupdates(rateinfo)))
        else
            # No more events
            break
        end
    end

    # Make sure all amounts are exactly computed at end
    final_time = isfinite(endtime) ? endtime : cur_time

    for (trap, cur_fill) ∈ enumerate(cur_amounts)
        if cur_fill.time < endtime
            if leakage_state.leaking[trap]
                # For leaking trap, use dynamic equilibrium volume (accounts for
                # drainage only during periods without inflow)
                drained_vol = compute_dynamic_equilibrium_volume(trap, final_time, leakage_state)
                # For leaking traps, drained_vol should never be nothing
                final_vol = isnothing(drained_vol) ? cur_fill.amount : drained_vol
                cur_amounts[trap] = FilledAmount(final_vol, final_time)
            else
                cur_amounts[trap] = FilledAmount(
                    SurfaceWaterIntegratedModeling._compute_exact_fill(
                        rateinfo, cur_amounts, trap,
                        filled_traps, tstruct, endtime,
                        z_vol_tables, false
                    ),
                    final_time
                )
            end
        end
    end
end


"""
    get_effective_leakage_cap(leakage_state, trap_id) -> Float64

Get the effective volume cap for a leaking trap.

Returns the leakage_volume for traps that are leaking (volume at which CO2 height
reaches the leakage threshold). Returns Inf for non-leaking traps.

NOTE: Only traps that actually started leaking are marked as leaking. Ancestor traps
are NOT marked as leaking - they continue accumulating normally and spill to the
leaking trap, which then passes CO2 out of domain. This allows multiple injection
sites to work correctly when trap hierarchies share common ancestors.
"""
function get_effective_leakage_cap(leakage_state::LeakageState, trap_id::Int)::Float64
    if !leakage_state.leaking[trap_id]
        # Not leaking - no cap (return Inf)
        return Inf
    end

    # This trap is leaking - cap at its leakage volume
    return leakage_state.leakage_volume[trap_id]
end


"""
    _fix_leaking_trap_edges!(sgraph, graph_updates, leakage_state)

Restore edge=0 for all leaking traps after `update_spillgraph!`.

SWIM's `update_spillgraph!` contains a sibling cycle check that redirects all
sibling traps to their parent when all become "full". This can overwrite edge=0
for previously leaking traps, breaking flow routing. This helper restores the
correct edges and ensures `graph_updates` reflects the fix so `_update_flow!`
processes these traps.
"""
function _fix_leaking_trap_edges!(sgraph, graph_updates, leakage_state::LeakageState)
    for trap in eachindex(leakage_state.leaking)
        if leakage_state.leaking[trap] && get(sgraph.edges, trap, 0) != 0
            old_target = sgraph.edges[trap]
            sgraph.edges[trap] = 0
            # Also fix in graph_updates so _update_flow! processes this trap.
            # graph_updates values are (old_target, new_target) tuples.
            found = false
            for i in eachindex(graph_updates)
                if graph_updates[i].index == trap
                    prev_old = graph_updates[i].value[1]
                    graph_updates[i] = IncrementalUpdate(trap, (prev_old, 0))
                    found = true
                    break
                end
            end
            if !found
                push!(graph_updates, IncrementalUpdate(trap, (old_target, 0)))
            end
        end
    end
end
