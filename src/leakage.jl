using SurfaceWaterIntegratedModeling
import Graphs
import Interpolations
using DifferentialEquations: solve, ODEProblem, VectorContinuousCallback, terminate!

export get_leakage_root_trap
export fill_sequence_until_leak

function get_leakage_root_trap(
    tstruct::TrapStructure{<:Real},
    injection_location::Int,
    leakage_height::Float64
)::Union{Int,Nothing}
    N = numtraps(tstruct)

    spillpoint_heights_over_trap_bottom = Vector{Float64}(undef, N)
    for i = 1:N
        sp = tstruct.spillpoints[i]
        trap_bottom = minimum(tstruct.topography[tstruct.footprints[i]])
        spillpoint_heights_over_trap_bottom[i] = sp.elevation - trap_bottom
    end

    # identify the trap containing the injection location and that has a
    # spillpoint height over trap bottom that is at least leakage_height
    traps_higher_than_leakage = findall(
        x -> x > leakage_height,
        spillpoint_heights_over_trap_bottom
    )

    # Iterate over the traps higher that leakage, and find the one with the smallest spillpoint_heights_over_trap_bottom that contains the injection location
    candidate_traps = Int[]
    for trap in traps_higher_than_leakage
        if injection_location in tstruct.footprints[trap]
            push!(candidate_traps, trap)
        end

    end
    if !isempty(candidate_traps)
        # return the trap with the smallest spillpoint_heights_over_trap_bottom
        sorted_candidates = sort(candidate_traps, by=x -> spillpoint_heights_over_trap_bottom[x])
        return sorted_candidates[1]
    end

    return nothing
end

function get_seq_index_when_children_are_filled(
    tstruct::TrapStructure{<:Real},
    seq::Vector{SpillEvent},
    trap_index::Int
)::Union{Int,Nothing}
    children_of = x -> Graphs.inneighbors(tstruct.agglomerations, x)
    children_of_trap = children_of(trap_index)
    filled_idx = findfirst(i -> all(filled_at(seq, i)[children_of_trap]), 1:length(seq))
    return filled_idx
end

"""
    compute_complete_spillgraph(tstruct, full_traps, root_trap=nothing)

Compute the complete spillgraph, given a trapping structure and a (sub)set of
filled traps.

Returns an object of type [`SpillGraph`](@ref).  

The `full_traps` vector must be *consistent*, i.e. no filled trap can have
non-filled subtraps.

# Arguments
- `tstruct::TrapStructure{T<:Real}`: trapping structure
- `full_traps::Vector{Bool}`: a vector with one entry per trap, specifying if that
                              trap has been completely filled or not.
- `root_trap::Union{Int,Nothing}`: optional maximum parent trap. If specified,
                              traps will not spill beyond this parent trap.

See also [`TrapStructure`](@ref), [`update_spillgraph!`](@ref).
"""
function compute_complete_spillgraph_up_to_root(tstruct::TrapStructure{T},
    full_traps::Vector{Bool},
    root_trap::Union{Int,Nothing}=nothing) where T<:Real

    num_traps = numtraps(tstruct)
    @assert num_traps == length(full_traps)

    # defining two simple shorthands
    children_of = x -> Graphs.inneighbors(tstruct.agglomerations, x)
    parent_of = x -> Graphs.outneighbors(tstruct.agglomerations, x)

    # identify the traps that are filled 
    @assert SurfaceWaterIntegratedModeling._valid_trap_status(full_traps, tstruct)

    # identify traps that are active 
    active_traps = [full_traps[x] | all(full_traps[children_of(x)]) for x in 1:num_traps]

    # determine where the filled traps spill
    full_trap_ixs = findall(full_traps)

    function target(x)
        parent = parent_of(x)
        ds_reg = tstruct.spillpoints[x].downstream_region
        if !isempty(parent) && active_traps[parent[1]]
            # Check if we should respect root_trap constraint
            if !isnothing(root_trap) && parent[1] == root_trap
                # Don't spill beyond root_trap, treat as out-of-domain
                return num_traps + 1
            end
            return parent[1]
        end
        return ds_reg > 0 ? ds_reg : num_traps + 1 # latter value flags out-of-domain
    end

    spillgraph = SpillGraph()
    for x in full_trap_ixs
        SurfaceWaterIntegratedModeling._set_outedge!(spillgraph, x, target(x))
    end

    return spillgraph
end

"""
    update_spillgraph!(spillgraph, fill_changes, tstruct, root_trap=nothing)

Update an existing spillgraph after there have been changes in which traps are
filled or not.

The argument `spillgraph` will be modified.

# Arguments
- `spillgraph::SpillGraph`: the spillgraph to be updated
- `fill_changes::Vector{IncrementalUpdate{Bool}}`: changes in the fill status of 
       one or more traps, presented as a vector of incremental updates
- `tstruct::TrapStructure{T<:Real}`: the [`TrapStructure`](@ref) representing the
       trapping structure of the terrain the spill graph is associated with
- `root_trap::Union{Int,Nothing}`: optional root trap. If specified,
                              traps will not spill beyond this root trap.

See also [`IncrementalUpdate`](@ref), [`SpillGraph`](@ref) and 
[`compute_complete_spillgraph`](@ref).
"""
function update_spillgraph_with_root!(spillgraph::SpillGraph,
    fill_changes::Vector{IncrementalUpdate{Bool}},
    tstruct::TrapStructure{T},
    root_trap::Union{Int,Nothing}=nothing) where T<:Real

    num_traps = numtraps(tstruct)
    changes = Vector{Tuple{Int,Int,Int}}() # source, old target, new target

    for fc in fill_changes
        if fc.value # filled trap
            dsreg = tstruct.spillpoints[fc.index].downstream_region
            # num_traps+1 represents 'out of domain'
            target = dsreg > 0 ? dsreg : num_traps + 1
            push!(changes, SurfaceWaterIntegratedModeling._set_outedge!(spillgraph, fc.index, target))
        else # un-filled trap
            push!(changes, SurfaceWaterIntegratedModeling._remove_outedge!(spillgraph, fc.index))
        end
    end

    # Now that all updates have been initially applied, check for sibling cycles
    # that must be updated.  (We do this in a second pass, since it is not
    # guaranteed that the `fill_changes` have been presented in chronological order).
    for fc in fill_changes

        # identify siblings (if any)
        siblings, parent = SurfaceWaterIntegratedModeling._siblings_of(tstruct, fc.index)

        isempty(siblings) && continue # if the trap has no sibling, nothing to do

        if fc.value
            # This trap just filled.  Check if it completed a 'cycle' so that
            # flow should be redirected to its parent
            if all([SurfaceWaterIntegratedModeling._isfull(x, spillgraph) for x in siblings])
                # all siblings are full.  Ensure they are all spilling to parent
                for s in siblings
                    # Check if we should respect root_trap constraint
                    target_parent = (!isnothing(root_trap) && parent == root_trap) ?
                                    num_traps + 1 : parent
                    push!(changes, SurfaceWaterIntegratedModeling._set_outedge!(spillgraph, s, target_parent))
                end
            end
        else
            # This trap just stopped being filled.  Redirect any flow from other
            # siblings that still flow into parent.
            for s in siblings
                if s ∈ keys(spillgraph.edges)
                    @assert spillgraph.edges[s] == parent
                    dsreg = tstruct.spillpoints[s].downstream_region
                    target = dsreg > 0 ? dsreg : num_traps + 1
                    push!(changes, SurfaceWaterIntegratedModeling._set_outedge!(spillgraph, s, target))
                end
            end
        end
    end

    # Make the list of incremental updates.  If there are multiple changes
    # registered for a trap, merge them to ensure the old and new targets are
    # consistent with the old and new graph
    tmp = Dict{Int,Tuple{Int,Int}}()
    for ch in changes
        existing = get(tmp, ch[1], nothing)
        if isnothing(existing)
            tmp[ch[1]] = (ch[2], ch[3]) # from/to
        else
            tmp[ch[1]] = (existing[1], ch[3]) # re-use the old from value
        end
    end
    return [IncrementalUpdate(k, v) for (k, v) ∈ tmp]
end



function _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
    filled_traps, cur_amounts, z_vol_tables,
    tstruct, infiltration, endtime, time_slack,
    root_trap, verbose)
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

        # Check if root_trap has been filled - if so, terminate sequence
        if !isnothing(root_trap) && filled_traps[root_trap]
            break
        end
        # given changes in fill state, update spill graph
        graph_updates = update_spillgraph_with_root!(sgraph, fill_updates, tstruct, root_trap)

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

function fill_sequence_until_leak(tstruct::TrapStructure{<:Real},
    weather_events::Vector{WeatherEvent};
    time_slack::Real=0.0,
    infiltration::Union{Matrix{<:Real},Nothing}=nothing,
    root_trap::Union{Int,Nothing}=nothing,
    leakage_height::Union{Float64,Nothing}=nothing,
    verbose::Bool=false)::Vector{SpillEvent}
    @assert !isempty(weather_events)


    num_traps = numtraps(tstruct)
    (num_traps == 0) && return # if the terrain has no traps, there is nothing to do

    # initialize infiltration map from user input
    infiltration =
        (typeof(infiltration) == Nothing) ? zeros(size(tstruct.topography)) :
        (typeof(infiltration) <: Real) ? ones(size(tstruct.topography)) * infiltration :
        infiltration
    # compute tables to support computation of trap water volume as function of
    # water level
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), num_traps)
    sgraph = compute_complete_spillgraph_up_to_root(tstruct, filled_traps, root_trap)

    # start with empty sequence
    seq = Vector{SpillEvent}()

    # compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp

        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # compute inflow/runoff/infiltration rates corresponding to the fill
        # graph and new rain rate
        rateinfo = compute_flow(sgraph, we.rain_rate, infiltration, tstruct, verbose)

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
        _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
            filled_traps, cur_amounts, z_vol_tables,
            tstruct, infiltration, end_time, time_slack,
            root_trap, verbose)

        # if !isnothing(root_trap) && !isnothing(leakage_height) && !filled_traps[root_trap]
        #     # Now we have filled the seq up so only the root_trap remains to be filled
        #     # We now want to fill the root_trap until it reaches the leakage_height.
        #     # Register this a last event in the sequence
        # end

        if !isnothing(root_trap) && !isnothing(leakage_height) && !filled_traps[root_trap]
            # Now we have filled the seq up so only the root_trap remains to be filled
            # We now want to fill the root_trap until it reaches the leakage_height.
            # Register this a last event in the sequence
            # Compute target z and corresponding volume using z_vol_tables
            zv = z_vol_tables[root_trap]
            zvals, vvals = zv[1], zv[2]
            footprint = tstruct.footprints[root_trap]
            trap_bottom = tstruct.topography[footprint]
            zmin = minimum(trap_bottom)
            ztarget = zmin + leakage_height
            z2v = length(zvals) == 1 ?
                  (z -> 0.0) :
                  Interpolations.linear_interpolation(zvals, vvals,
                extrapolation_bc=Interpolations.Line())
            tvolume = (root_trap > numregions(tstruct)) ?
                      tstruct.trapvolumes[root_trap] - tstruct.subvolumes[root_trap] :
                      tstruct.trapvolumes[root_trap]
            vtarget = clamp(z2v(ztarget), 0.0, tvolume)

            cur_amt = cur_amounts[root_trap]
            if cur_amt.amount < vtarget
                vol_end, dt = fill_trap_until_volume(root_trap, rateinfo, cur_amt,
                    Inf, tstruct, z_vol_tables, vtarget, use_saved=false)
                if !isnothing(dt)
                    # # WARNING: The cur_amt.time is Inf. We should probably set it to the dt
                    # global_time = cur_amt.time + dt

                    time_at_prev_fill = seq[length(seq)].timestamp
                    println("time at prev fill: ", time_at_prev_fill)
                    println("dt: ", dt)
                    global_time = time_at_prev_fill + dt # FIX: Temporary hack

                    # update amounts and register one final event
                    cur_amounts[root_trap] = FilledAmount(vol_end, global_time)
                    amount_updates = [IncrementalUpdate(root_trap, FilledAmount(vol_end, global_time))]
                    fill_updates = Vector{IncrementalUpdate{Bool}}()
                    push!(seq, SpillEvent(global_time, amount_updates, fill_updates,
                        getinflowupdates(rateinfo), nothing, getrunoffupdates(rateinfo)))
                end
            end
        end

    end

    return seq
end

function fill_trap_until_volume(trap, rateinfo, cur_amount, endtime, tstruct, z_vol_tables, vtarget;
    use_saved=false)

    footprint = tstruct.footprints[trap]
    trap_bottom = tstruct.topography[footprint]
    tvolume = tstruct.trapvolumes[trap] - tstruct.subvolumes[trap]

    Smin = use_saved ? getsavedsmin(rateinfo, trap) : getsmin(rateinfo, trap)
    Smax = use_saved ? getsavedsmax(rateinfo, trap) : getsmax(rateinfo, trap)
    inflow = use_saved ? getsavedinflow(rateinfo, trap) : getinflow(rateinfo, trap)

    if (trap > numregions(tstruct))
        children = subtrapsof(tstruct, trap)
        # WARNING: Here we use the height to the top of the spillpoint, not the bottom
        # trap_bottom = max.(trap_bottom, tstruct.spillpoints[children[1]].elevation)
        trap_bottom = trap_bottom
    else
        tvolume = tstruct.trapvolumes[trap]
    end

    # clamp vtarget to feasible range
    vtarget = clamp(vtarget, 0.0, tvolume)

    if Smax == Smin
        # WARNING: For us the inflow is always zero (as the filling is from a well)
        # so this case is always hit
        # accum_rate = inflow
        accum_rate = 1.0 # FIX: Temorary
        if accum_rate == 0.0
            return (cur_amount.amount, nothing)
        end
        # time to reach vtarget (not necessarily full/empty)
        if accum_rate > 0.0
            if vtarget <= cur_amount.amount
                return (cur_amount.amount, cur_amount.time) # already at/above target
            end
            dt_target = (vtarget - cur_amount.amount) / accum_rate
            t_target = cur_amount.time + dt_target
            reached = (t_target <= endtime)
            return reached ? (vtarget, dt_target) :
                   (cur_amount.amount + (endtime - cur_amount.time) * accum_rate, nothing)
        else
            # decreasing volume: if current volume already <= vtarget then already reached
            if cur_amount.amount <= vtarget
                return (cur_amount.amount, cur_amount.time)
            end
            dt_target = (cur_amount.amount - vtarget) / abs(accum_rate)
            t_target = cur_amount.time + dt_target
            reached = (t_target <= endtime)
            return reached ? (vtarget, dt_target) :
                   (max(0.0, cur_amount.amount + (endtime - cur_amount.time) * accum_rate), nothing)
        end
    end

    # use ODE integration otherwise
    fprint_infil =
        use_saved ? [-min(getsavedrunoff(rateinfo, i), 0.0) for i in footprint] :
        -min.(getrunoff(rateinfo, footprint), 0.0)

    infilfun = ixs -> sum(fprint_infil[ixs]) - Smin
    v0 = [cur_amount.amount]
    dvdt = SurfaceWaterIntegratedModeling._setup_dvdt(trap_bottom, tvolume, infilfun, inflow,
        tstruct.spillpoints[trap], z_vol_tables[trap])
    dv0 = [0.0]
    dvdt(dv0, v0, 0, 0)

    (v0[1] == 0.0 && dv0[1] <= 0) && return (0.0, cur_amount.time)
    (v0[1] == tvolume && dv0[1] >= 0) && return (tvolume, cur_amount.time)
    if v0[1] >= vtarget && dv0[1] <= 0
        return (v0[1], cur_amount.time)
    end

    function condition(out, v, t, integrator)
        out[1] = v[1] - vtarget    # reach target volume
        out[2] = tvolume - v[1]    # full
        out[3] = dv0[1] * begin
            deriv = [0.0]
            dvdt(deriv, v, 0, t)
            deriv[1]
        end
    end

    condition_reached = [0]
    function affect!(integrator, ix)
        condition_reached[] = ix
        terminate!(integrator)
    end
    cb = VectorContinuousCallback(condition, affect!, 3)
    dt = endtime - cur_amount.time
    sol = solve(ODEProblem(dvdt, v0, [0, dt]), callback=cb, abstol=1e-6, reltol=1e-4)

    return (sol.u[end][1],
        (condition_reached[1] ∈ [1, 2]) ? sol.t[end] : nothing)
end

