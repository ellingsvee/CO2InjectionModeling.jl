using SurfaceWaterIntegratedModeling
import Interpolations
using DifferentialEquations: solve, ODEProblem, VectorContinuousCallback, terminate!
export fill_layer, InjectionEvent


function fill_layer(tstruct::TrapStructure{<:Real},
            domain::Domain3D,
            reservoir_properties::ReservoirProperties,
            weather_events::Vector{WeatherEvent};
            time_slack::Float64=0.0, # NOT USED. Included for legacy
            infiltration::Union{Matrix{<:Real}, Nothing} = nothing,
            verbose::Bool=false)
    @assert !isempty(weather_events)

    num_traps = SurfaceWaterIntegratedModeling.numtraps(tstruct)
    (num_traps == 0) && return # if the terrain has no traps, there is nothing to do

    # initialize infiltration map from user input
    infiltration =
        (typeof(infiltration) == Nothing) ? zeros(size(tstruct.topography)) :
        (typeof(infiltration) <: Real)  ? ones(size(tstruct.topography)) * infiltration :
                                          infiltration
    # compute tables to support computation of trap water volume as function of
    # water level
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)


    # set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), num_traps)
    sgraph = SurfaceWaterIntegratedModeling.compute_complete_spillgraph(tstruct, filled_traps)

    # Initialize leakage state
    leakage_threshold = compute_leakage_height(reservoir_properties)
    leakage_state = LeakageState(
        falses(num_traps),              # leaked_traps
        Vector{LeakageEvent}(),         # leakage_events
        leakage_threshold,              # leakage_threshold
        nothing                         # first_leakage_time
    )

    # start with empty sequence
    seq = Vector{SurfaceWaterIntegratedModeling.SpillEvent}()
    snapshots = Vector{SimulationLayerSnapshot}()

    # Initialize source tracker with the first weather event's rain_rate
    source_tracker = SourceTracker(tstruct, weather_events[1].rain_rate)

    # compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp

        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # Make a mutable copy of rain_rate for leakage handling
        rain_rate = copy(we.rain_rate)

        # Apply leakage state to zero out injection in regions that have already leaked
        # This ensures that once injection is stopped due to leakage, it stays stopped
        # across weather event boundaries
        _apply_leakage_to_rain_rate!(rain_rate, leakage_state, tstruct)

        # Update source tracker with new injection sources from this weather event
        # This is important for multi-layer simulations where leakage from upper layers
        # creates new injection points in lower layers
        # NOTE: Use the modified rain_rate (after applying leakage) to avoid tracking
        # sources that have been shut off due to leakage
        update_injection_sources!(source_tracker, tstruct, rain_rate)

        # compute inflow/runoff/infiltration rates corresponding to the fill
        # graph and new rain rate
        # IMPORTANT: Use the modified rain_rate (after applying leakage) so that
        # injection in leaked regions is not included in flow computation
        rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, rain_rate, infiltration, tstruct, verbose)

        # compute initial time estimates for when a trap become filled, or split
        # into subtraps
        changetimeest = SurfaceWaterIntegratedModeling._set_initial_changetime_estimates(rateinfo, cur_amounts,
                                                          cur_time, filled_traps,
                                                          tstruct)

        # register the start of this weather event as a new, fully computed, spill event
        push!(seq, SurfaceWaterIntegratedModeling.SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
                              copy(rateinfo.trap_inflow), copy(rain_rate),
                              copy(rateinfo.runoff)))

        # Will add new events to `seq`.  `sgraph`, `rateinfo`, `changetimeest`,
        # `filled_traps` and `cur_amounts` are also modified in the process
        _fill_sequence_for_weather_event_with_leakage!(seq, sgraph, rateinfo, changetimeest,
                                          filled_traps, cur_amounts, z_vol_tables,
                                          tstruct, infiltration, end_time, time_slack,
                                          leakage_state, rain_rate, source_tracker,
                                          reservoir_properties, verbose)
    end

    return seq, leakage_state
end

"""
Helper function to identify all parent traps containing a given trap.
For regions (trap_id <= numregions), returns all containing parent traps.
For parent traps, returns itself and any higher-level parents.
"""
function _identify_all_parent_traps(trap_id::Int, tstruct)
    parent = parentof(tstruct, trap_id)
    if isnothing(parent)
        return Int[]
    else
        return [parent; _identify_all_parent_traps(parent, tstruct)...]
    end
end

"""
Find all filled parent traps that are connected to a leaking trap.

These are the traps whose CO2 will drain through the leaking trap.
Only filled parent traps contribute their volume to the leakage.
"""
function _find_connected_filled_parents(trap_id::Int, filled_traps::Vector{Bool}, tstruct)
    connected = Int[trap_id]  # Start with the leaking trap itself

    # Traverse up the parent hierarchy
    current = trap_id
    while true
        parent = parentof(tstruct, current)
        if isnothing(parent)
            break
        end

        # Only include if the parent is filled
        if filled_traps[parent]
            push!(connected, parent)
            current = parent
        else
            # If parent is not filled, stop traversing
            break
        end
    end

    return connected
end


function _fill_sequence_for_weather_event_with_leakage!(seq, sgraph, rateinfo, changetimeest,
                                           filled_traps, cur_amounts, z_vol_tables,
                                           tstruct, infiltration, endtime, time_slack,
                                           leakage_state::LeakageState, rain_rate::Matrix{Float64},
                                           source_tracker::SourceTracker,
                                           reservoir_properties::ReservoirProperties, verbose)
    cur_time = cur_amounts[1].time

    fill_updates = Vector{IncrementalUpdate{Bool}}()
    graph_updates = Vector{IncrementalUpdate{Int}}()

    count = 0
    while cur_time < endtime
        verbose && (mod(count+=1, 10) == 0) && println("Fill sequence iteration: ", count)

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

        # Update source tracking for newly filled traps
        # This must happen AFTER the spill graph is updated so we know where each trap spills to
        for u in fill_updates
            if u.value  # trap just filled (not emptied)
                update_sources_on_trap_fill!(source_tracker, u.index, tstruct, sgraph)
            end
        end

        # given the updates ot the spill graph, update flow information in `rateinfo`
        SurfaceWaterIntegratedModeling.setsavepoint!(rateinfo)
        SurfaceWaterIntegratedModeling._update_flow!(rateinfo, graph_updates, tstruct, sgraph)

        # update water amount in traps whose inflow rate is about to change, or
        # that just filled
        amount_updates = SurfaceWaterIntegratedModeling._update_affected_amounts(rateinfo, cur_amounts, filled_traps,
                                                  tstruct, z_vol_tables, cur_time)
        append!(amount_updates,
                [IncrementalUpdate(tix, SurfaceWaterIntegratedModeling.FilledAmount(tstruct.trapvolumes[tix] -
                    tstruct.subvolumes[tix], cur_time))
                 for tix in [u.index for u in fill_updates]])

        # integrate the changes into the continously updated `cur_amounts` vector
        SurfaceWaterIntegratedModeling._apply_updates!(cur_amounts, amount_updates)

        # add current state to result
        push!(seq, SurfaceWaterIntegratedModeling.SpillEvent(cur_time, amount_updates, fill_updates,
                              SurfaceWaterIntegratedModeling.getinflowupdates(rateinfo), nothing,
                              SurfaceWaterIntegratedModeling.getrunoffupdates(rateinfo)))

        # CHECK FOR LEAKAGE after each status change
        leakage_occurred = _check_and_trigger_leakage!(leakage_state, cur_amounts,
                                                       rain_rate, z_vol_tables,
                                                       tstruct, cur_time, source_tracker,
                                                       filled_traps, reservoir_properties, verbose, rateinfo)

        # If leakage occurred, we need to recompute flow rates with modified injection
        if leakage_occurred
            # Recompute flow with updated rain_rate (zeroed out in leaked regions)
            SurfaceWaterIntegratedModeling.setsavepoint!(rateinfo)
            rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, rain_rate,
                                                                   infiltration, tstruct, verbose)

            # Add a new SpillEvent to record the state after leakage
            # This is critical for correct interpolation in trap_states_at_timepoints
            # Without this event, interpolation will use stale flow rates from before leakage
            push!(seq, SurfaceWaterIntegratedModeling.SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
                                  copy(rateinfo.trap_inflow), copy(rain_rate),
                                  copy(rateinfo.runoff)))

            # Recompute time estimates with new flow rates
            changetimeest = SurfaceWaterIntegratedModeling._set_initial_changetime_estimates(
                rateinfo, cur_amounts, cur_time, filled_traps, tstruct)
        end
    end

    # make sure all amounts are exactly computed at end
    for (trap, cur_fill) ∈ enumerate(cur_amounts)
        if cur_fill.time < endtime
            cur_amounts[trap] =
                SurfaceWaterIntegratedModeling.FilledAmount(SurfaceWaterIntegratedModeling._compute_exact_fill(rateinfo, cur_amounts, trap,
                                                 filled_traps, tstruct, endtime,
                                                 z_vol_tables, false),
                             min(cur_time, endtime))
        end
    end
end