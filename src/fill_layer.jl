using SurfaceWaterIntegratedModeling
import Interpolations

export fill_sequence_with_leakage

"""
Modified version of the SWIM fill sequence algorithm that accounts for CO2 leakage.
"""
function fill_sequence_with_leakage(tstruct::TrapStructure{<:Real},
    weather_events::Vector{WeatherEvent};
    time_slack::Real=0.0,
    verbose::Bool=false)::Vector{SpillEvent}
    @assert !isempty(weather_events)

    num_traps = numtraps(tstruct)
    (num_traps == 0) && return # if the terrain has no traps, there is nothing to do

    # compute tables to support computation of trap water volume as function of
    # water level
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), num_traps)
    sgraph = SurfaceWaterIntegratedModeling.compute_complete_spillgraph(tstruct, filled_traps)

    # start with empty sequence
    seq = Vector{SpillEvent}()

    # compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp

        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # compute inflow/runoff rates corresponding to the fill
        # graph and new rain rate
        rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, we.rain_rate, zeros(size(tstruct.topography)), tstruct, verbose)

        # compute initial time estimates for when a trap become filled, or split
        # into subtraps
        changetimeest = SurfaceWaterIntegratedModeling._set_initial_changetime_estimates(rateinfo, cur_amounts,
            cur_time, filled_traps,
            tstruct)

        # register the start of this weather event as a new, fully computed, spill event
        push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
            copy(rateinfo.trap_inflow), copy(we.rain_rate),
            copy(rateinfo.runoff)))

        # Will add new events to `seq`.  `sgraph`, `rateinfo`, `changetimeest`,
        # `filled_traps` and `cur_amounts` are also modified in the process
        _fill_sequence_with_leakage_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
            filled_traps, cur_amounts, z_vol_tables,
            tstruct, end_time, time_slack,
            verbose)
    end

    return seq
end

function _fill_sequence_with_leakage_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
    filled_traps, cur_amounts, z_vol_tables,
    tstruct, endtime, time_slack,
    verbose)
    cur_time = cur_amounts[1].time

    fill_updates = Vector{IncrementalUpdate{Bool}}()
    graph_updates = Vector{IncrementalUpdate{Int}}()

    count = 0
    while cur_time < endtime
        verbose && (mod(count += 1, 10) == 0) && println("Fill sequence iteration: ", count)

        cur_time, fill_updates =
            SurfaceWaterIntegratedModeling._identify_next_status_change!(changetimeest, cur_amounts, rateinfo,
                filled_traps, tstruct, z_vol_tables,
                cur_time, endtime)

        (cur_time > endtime || isempty(fill_updates)) && break # do not register
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
    end

    # make sure all amounts are exactly computed at end
    for (trap, cur_fill) ∈ enumerate(cur_amounts)
        if cur_fill.time < endtime
            cur_amounts[trap] =
                FilledAmount(SurfaceWaterIntegratedModeling._compute_exact_fill(rateinfo, cur_amounts, trap,
                        filled_traps, tstruct, endtime,
                        z_vol_tables, false),
                    min(cur_time, endtime))
        end
    end
end