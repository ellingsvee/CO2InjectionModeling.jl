using SurfaceWaterIntegratedModeling
import Interpolations

export fill_sequence_with_leakage

"""
Modified version of the SWIM fill sequence algorithm that accounts for CO2 leakage.
"""
function fill_sequence_with_leakage(tstruct::TrapStructure{<:Real},
    reservoir_properties::ReservoirProperties,
    weather_events::Vector{WeatherEvent};
    time_slack::Real=0.0,
    verbose::Bool=false)
    @assert !isempty(weather_events)

    num_traps = numtraps(tstruct)
    # (num_traps == 0) && return # if the terrain has no traps, there is nothing to do
    if num_traps == 0
        empty_leakage_state = create_empty_leakage_state()
        return Vector{SpillEvent}(), empty_leakage_state
    end

    # compute tables to support computation of trap water volume as function of
    # water level
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), num_traps)
    sgraph = SurfaceWaterIntegratedModeling.compute_complete_spillgraph(tstruct, filled_traps)

    # start with empty sequence
    seq = Vector{SpillEvent}()

    # initialize leakage state
    leakage_state = initialize_leakage_state(
        tstruct, z_vol_tables, reservoir_properties,
        reservoir_properties.sand_residual_co2_saturation,
        reservoir_properties.residual_leakage_time
    )

    # compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp

        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # compute inflow/runoff rates corresponding to the fill
        # graph and new rain rate
        rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, we.rain_rate, zeros(size(tstruct.topography)), tstruct, verbose)

        # compute initial time estimates for when a trap become filled, or split into subtraps
        changetimeest = SurfaceWaterIntegratedModeling._set_initial_changetime_estimates(rateinfo, cur_amounts,
            cur_time, filled_traps, tstruct)

        # initial estimate of leakage times
        leakage_time_est = get_initial_leakage_time_estimates(cur_amounts, cur_time, rateinfo, leakage_state, filled_traps, tstruct)

        # register the start of this weather event as a new, fully computed, spill event
        push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
            copy(rateinfo.trap_inflow), copy(we.rain_rate),
            copy(rateinfo.runoff)))

        #  will add new events to `seq`.  `sgraph`, `rateinfo`, `changetimeest`,
        # `filled_traps` and `cur_amounts` are also modified in the process
        _fill_sequence_with_leakage_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
            filled_traps, cur_amounts, z_vol_tables,
            tstruct, end_time, time_slack, leakage_state, leakage_time_est,
            verbose)
    end

    return seq, leakage_state
end

function _fill_sequence_with_leakage_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
    filled_traps, cur_amounts, z_vol_tables,
    tstruct, endtime, time_slack, leakage_state, leakage_time_est,
    verbose)

    cur_time = cur_amounts[1].time
    num_traps = numtraps(tstruct)

    fill_updates = Vector{IncrementalUpdate{Bool}}()
    graph_updates = Vector{IncrementalUpdate{Int}}()

    count = 0
    while cur_time < endtime
        verbose && (mod(count += 1, 10) == 0) && println("Fill sequence iteration: ", count)

        next_fill_time, fill_updates =
            SurfaceWaterIntegratedModeling._identify_next_status_change!(changetimeest, cur_amounts, rateinfo,
                filled_traps, tstruct, z_vol_tables,
                cur_time, endtime)

        # Find next leakage event
        next_leak_time, leak_trap = find_next_leakage_event(leakage_time_est)

        if next_leak_time < next_fill_time && next_leak_time <= endtime && leak_trap > 0
            cur_time = next_leak_time

            verbose && println("Leakage event at time $(cur_time) in trap $(leak_trap)")

            # Mark trap as leaking
            leakage_state.leaking[leak_trap] = true

            # Only mark as draining if trap has actual volume (not pass-through with vol 0)
            leakage_vol_check = leakage_state.leakage_volume[leak_trap]
            if leakage_vol_check > 0
                leakage_state.draining[leak_trap] = true
                leakage_state.leakage_start_time[leak_trap] = cur_time
            end

            # Record the leakage for upstream layer
            leakage_location = find_leakage_location(leak_trap, tstruct)
            push!(leakage_state.leakage_records, LeakageRecord(
                cur_time,
                leak_trap,
                leakage_location
            ))

            # Record the initial volume at the time leakage started (for residual drainage)
            leakage_vol = leakage_state.leakage_volume[leak_trap]
            leakage_state.initial_volume_at_leak[leak_trap] = leakage_vol

            # Update desendants of the leaking trap
            descendants = get_all_descendants(tstruct, leak_trap)
            for desc_id in descendants
                @assert filled_traps[desc_id] == true "Descendant trap $(desc_id) is not filled at time of leakage"
                @assert leakage_state.leaking[desc_id] == false "Descendant trap $(desc_id) is already leaking"
                @assert leakage_state.draining[desc_id] == false "Descendant trap $(desc_id) is already draining"

                leakage_state.leakage_start_time[desc_id] = cur_time
                leakage_state.draining[desc_id] = true

                desc_vol = cur_amounts[desc_id].amount
                leakage_state.initial_volume_at_leak[desc_id] = desc_vol

            end



            # Cap the trap amount at leakage volume
            cur_amounts[leak_trap] = FilledAmount(leakage_vol, cur_time)

            # Mark the trap as filled to prevent further fill events
            filled_traps[leak_trap] = true
            leakage_state.leaking[leak_trap] = true

            # Create fill update for the leaking trap only
            leak_fill_updates = [IncrementalUpdate{Bool}(leak_trap, true)]

            # Update spillgraph for the leaking trap
            graph_updates = SurfaceWaterIntegratedModeling.update_spillgraph!(sgraph, leak_fill_updates, tstruct)

            # TODO: Verify this logic...
            # Set ONLY the leaking trap's edge to 0 (out of domain = leakage)
            # Ancestors keep their normal edges - they spill to this trap, which then leaks
            sgraph.edges[leak_trap] = 0

            # Update flow information with the modified spillgraph
            setsavepoint!(rateinfo)
            SurfaceWaterIntegratedModeling._update_flow!(rateinfo, graph_updates, tstruct, sgraph)

            # Create amount update for this trap
            amount_updates = [IncrementalUpdate(leak_trap, FilledAmount(leakage_vol, cur_time))]

            # Update leakage time estimate to Inf for the leaking trap (already leaking)
            leakage_time_est[leak_trap] = LeakageTimeEstimate(leak_trap, Inf, Inf)


            # Update leakage time estimates for all traps whose inflow changed
            affected_traps = unique([u.index for u in getinflowupdates(rateinfo)])
            update_leakage_time_estimates!(
                leakage_time_est, affected_traps, cur_amounts, cur_time,
                rateinfo, leakage_state, filled_traps, tstruct
            )

            # Also update SWIM's changetime estimates for affected traps and the leaking trap
            all_traps_to_update = unique(vcat(affected_traps, [leak_trap]))
            for trap in all_traps_to_update
                changetimeest[trap] = SurfaceWaterIntegratedModeling._compute_changetime_estimate(
                    trap, cur_amounts, cur_time, rateinfo, filled_traps, tstruct
                )
            end

            # Record this as a spill event
            push!(seq, SpillEvent(cur_time, amount_updates, leak_fill_updates,
                getinflowupdates(rateinfo), nothing,
                getrunoffupdates(rateinfo)))

        elseif next_fill_time <= endtime && !isempty(fill_updates)
            cur_time = next_fill_time

            # more events
            for u in fill_updates
                filled_traps[u.index] = u.value
            end
            # given changes in fill state, update spill graph
            graph_updates = SurfaceWaterIntegratedModeling.update_spillgraph!(sgraph, fill_updates, tstruct)

            # given the updates ot the spill graph, update flow information in `rateinfo`
            setsavepoint!(rateinfo)
            SurfaceWaterIntegratedModeling._update_flow!(rateinfo, graph_updates, tstruct, sgraph)

            # update water amount in traps whose inflow rate is about to change, or
            # that just filled
            amount_updates = SurfaceWaterIntegratedModeling._update_affected_amounts(rateinfo, cur_amounts, filled_traps,
                tstruct, z_vol_tables, cur_time)

            append!(amount_updates,
                [IncrementalUpdate(tix, FilledAmount(tstruct.trapvolumes[tix] -
                                                     tstruct.subvolumes[tix], cur_time))
                 for tix in [u.index for u in fill_updates]])

            # integrate the changes into the continously updated `cur_amounts` vector
            SurfaceWaterIntegratedModeling._apply_updates!(cur_amounts, amount_updates)

            # add current state to result
            push!(seq, SpillEvent(cur_time, amount_updates, fill_updates,
                getinflowupdates(rateinfo), nothing,
                getrunoffupdates(rateinfo)))

        else
            break # no more events within this weather event
        end
    end


    # make sure all amounts are exactly computed at end
    for (trap, cur_fill) ∈ enumerate(cur_amounts)
        if cur_fill.time < endtime
            if leakage_state.leaking[trap]
                # Leaking trap - compute volume accounting for residual drainage
                # The stored volume decreases over time as CO2 drains out
                drained_vol = compute_volume_with_drainage(trap, endtime, leakage_state)

                # For leaking traps, drained_vol should never be nothing
                final_vol = isnothing(drained_vol) ? cur_fill.amount : drained_vol

                cur_amounts[trap] = FilledAmount(final_vol, endtime)
            else
                cur_amounts[trap] =
                    FilledAmount(SurfaceWaterIntegratedModeling._compute_exact_fill(rateinfo, cur_amounts, trap,
                            filled_traps, tstruct, endtime,
                            z_vol_tables, false),
                        endtime)
            end
        end
    end
end
