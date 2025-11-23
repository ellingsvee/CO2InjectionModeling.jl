using SurfaceWaterIntegratedModeling
import Interpolations
using DifferentialEquations: solve, ODEProblem, VectorContinuousCallback, terminate!
export fill_layer, InjectionEvent, LeakageState, LeakageEvent, compute_per_trap_leaked_volumes, get_total_leaked_volume

struct InjectionEvent
    timestamp::Float64
    injection_rate::Union{Matrix{Float64}, Float64} 
end

struct LeakageEvent
    timestamp::Float64           # Time when leakage was detected (at spill event)
    trap_id::Int                # Trap that leaked
    height_at_detection::Float64  # Height when leakage was detected (>= threshold)
    volume_in_trap::Float64     # Volume in trap at detection (swim units)
    source_regions::Set{Int}    # Regions whose injection contributed to this leakage
    connected_filled_traps::Vector{Int}  # Trap IDs that are filled and connected (will leak)
    total_leakable_volume::Float64       # Total volume that will leak (after residual trapping)
    residual_trapped_volume::Float64     # Volume that stays trapped (residual)
end

mutable struct LeakageState
    leaked_traps::Vector{Bool}              # Which traps have leaked
    leakage_events::Vector{LeakageEvent}    # Detailed leakage events
    leakage_threshold::Float64              # Height threshold (meters)
    first_leakage_time::Union{Float64, Nothing}  # Time of first leakage (Nothing if no leakage yet)
end

"""
Tracking structure that maps traps to their contributing injection sources.

This tracks the flow of CO2 from injection points through the trap hierarchy,
allowing us to determine which injection sources contributed to filling each trap.
"""
mutable struct SourceTracker
    # For each trap, store the set of source regions that contributed CO2 to it
    # Key: trap_id, Value: Set of region indices that contributed
    trap_sources::Dict{Int, Set{Int}}

    # Regions that have active injection (rain_rate > 0)
    injection_regions::Set{Int}
end

function SourceTracker(tstruct::TrapStructure, rain_rate::Matrix{<:Real})
    num_traps = SurfaceWaterIntegratedModeling.numtraps(tstruct)

    # Initialize trap_sources: each trap starts with empty source set
    trap_sources = Dict{Int, Set{Int}}()
    for t in 1:num_traps
        trap_sources[t] = Set{Int}()
    end

    # Find regions with active injection
    injection_regions = Set{Int}()
    for idx in CartesianIndices(rain_rate)
        if rain_rate[idx] > 0
            region = tstruct.regions[idx]
            if region > 0
                push!(injection_regions, region)
            end
        end
    end

    # Initialize: regions with injection are sources for themselves
    num_regions = SurfaceWaterIntegratedModeling.numregions(tstruct)
    for region in 1:num_regions
        if region in injection_regions
            push!(trap_sources[region], region)
        end
    end

    return SourceTracker(trap_sources, injection_regions)
end

"""
Update source tracking when a trap becomes filled.

When a trap fills, it inherits sources from:
1. Its own region (for lowest-level traps) - already set at initialization
2. Its subtraps (for parent traps) - CO2 from subtraps flows into parent
3. Upstream traps that spill into it - tracked when those traps fill and spill over

This function should be called when a trap fills, passing the spill graph
to determine where this trap's CO2 will flow next.
"""
function update_sources_on_trap_fill!(tracker::SourceTracker,
                                       trap_id::Int,
                                       tstruct::TrapStructure,
                                       sgraph::SurfaceWaterIntegratedModeling.SpillGraph)
    # Get sources from subtraps (for parent traps)
    children = SurfaceWaterIntegratedModeling.subtrapsof(tstruct, trap_id)
    for child in children
        union!(tracker.trap_sources[trap_id], tracker.trap_sources[child])
    end

    # Now propagate this trap's sources to its downstream trap
    downstream = get(sgraph.edges, trap_id, nothing)
    if !isnothing(downstream) && downstream <= SurfaceWaterIntegratedModeling.numtraps(tstruct)
        union!(tracker.trap_sources[downstream], tracker.trap_sources[trap_id])
    end
end

"""
Find all source regions that contributed CO2 to a leaked trap.

This traces back through the source tracking to find all regions
whose injection contributed to filling the leaked trap.
"""
function find_contributing_sources(tracker::SourceTracker, trap_id::Int)
    return tracker.trap_sources[trap_id]
end

"""
Find all regions that should have injection stopped when a trap leaks.

For each source region that contributed to the leak, we need to stop injection
in that region.
"""
function find_regions_to_stop(tracker::SourceTracker, trap_id::Int)
    contributing_sources = find_contributing_sources(tracker, trap_id)
    # The sources are region indices, so we return them directly
    return contributing_sources
end

"""
Update the source tracker when a new weather event starts with a different rain_rate.

This is needed for multi-layer simulations where the initial rain_rate may be zero
(no direct injection) but later weather events have non-zero rain_rate from leakage.
"""
function update_injection_sources!(tracker::SourceTracker, tstruct::TrapStructure, rain_rate::Matrix{<:Real})
    # Find regions with active injection in this weather event
    for idx in CartesianIndices(rain_rate)
        if rain_rate[idx] > 0
            region = tstruct.regions[idx]
            if region > 0
                # Add this region as an injection source
                push!(tracker.injection_regions, region)
                # Mark this region as a source for itself
                push!(tracker.trap_sources[region], region)
            end
        end
    end
end

function fill_layer(tstruct::TrapStructure{<:Real},
            domain::Domain3D,
            reservoir_properties::ReservoirProperties,
            weather_events::Vector{WeatherEvent};
            time_slack::Float64=0.0, # NOT USED. Included for legacy
            infiltration::Union{Matrix{<:Real}, Nothing} = nothing,
            verbose::Bool=false)
    @assert !isempty(weather_events)
    # # Convert injection events to weather events
    # weather_events = [WeatherEvent(ie.timestamp, physical_volume_to_swim_volume(ie.injection_rate, reservoir_properties, domain)) for ie in injection_events]


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

        # Update source tracker with new injection sources from this weather event
        # This is important for multi-layer simulations where leakage from upper layers
        # creates new injection points in lower layers
        update_injection_sources!(source_tracker, tstruct, we.rain_rate)

        # compute inflow/runoff/infiltration rates corresponding to the fill
        # graph and new rain rate
        rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, we.rain_rate, infiltration, tstruct, verbose)

        # compute initial time estimates for when a trap become filled, or split
        # into subtraps
        changetimeest = SurfaceWaterIntegratedModeling._set_initial_changetime_estimates(rateinfo, cur_amounts,
                                                          cur_time, filled_traps,
                                                          tstruct)

        # Make a mutable copy of rain_rate for leakage handling
        rain_rate = copy(we.rain_rate)

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



# """
# Compute per-trap leakage volumes for injection into the next layer.

# Returns a dictionary mapping trap_id -> total volume leaked from that trap (in swim units).
# This computes how much CO2 leaked through each trap after it exceeded the threshold.

# For each leaked trap, the leaked volume includes:
# - Injection directly into the leaked trap's footprint
# - Injection into any parent traps that contain this leaked trap

# This accounts for the fact that CO2 injected into a parent trap will drain down through
# its child traps and leak through whichever child reaches the threshold first.

# Parameters:
# - leakage_state: LeakageState from the simulation
# - injection_events: Vector of injection events used in the simulation
# - reservoir_properties: ReservoirProperties used in the simulation
# - domain: Domain3D used in the simulation
# - tstruct: TrapStructure for the layer
# - final_time: End time of the simulation (used to compute total leaked volume)

# Returns:
# - Dict{Int, Float64}: Maps trap_id to leaked volume in swim units
# """
# function compute_per_trap_leaked_volumes(leakage_state::LeakageState,
#                                          injection_events::Vector{InjectionEvent},
#                                          reservoir_properties::ReservoirProperties,
#                                          domain::Domain3D,
#                                          tstruct::TrapStructure,
#                                          final_time::Float64)
#     # If no leakage, return empty dictionary
#     if isnothing(leakage_state.first_leakage_time)
#         return Dict{Int, Float64}()
#     end

#     # Get trap-to-footprint mapping
#     nx, ny = size(tstruct.topography)

#     # Initialize result dictionary
#     leaked_volumes_per_trap = Dict{Int, Float64}()

#     # For each leaked trap, compute volume that leaked
#     for leakage_event in leakage_state.leakage_events
#         trap_id = leakage_event.trap_id
#         leak_time = leakage_event.timestamp

#         # Get all traps whose injection drains through this leaked trap
#         # This includes the trap itself AND all parent traps containing it
#         relevant_traps = _identify_all_parent_traps(trap_id, tstruct)
#         push!(relevant_traps, trap_id)

#         # Combine footprints from the leaked trap and all its parent traps
#         # The footprint of a trap represents all cells that drain into it
#         combined_footprint = Set{Int}()
#         for relevant_trap in relevant_traps
#             union!(combined_footprint, tstruct.footprints[relevant_trap])
#         end

#         # Sum injection into the combined footprint after leakage
#         total_leaked = 0.0

#         for (i, ie) in enumerate(injection_events)
#             start_time = ie.timestamp
#             end_time = (i < length(injection_events)) ? injection_events[i+1].timestamp : final_time

#             # If this injection period overlaps with post-leakage time
#             if end_time > leak_time
#                 # Compute how much time was spent injecting after leakage
#                 effective_start = max(start_time, leak_time)
#                 effective_duration = end_time - effective_start

#                 if effective_duration > 0
#                     # Convert injection rate to swim units
#                     rain_rate_swim = physical_volume_to_swim_volume(ie.injection_rate, reservoir_properties, domain)

#                     # Sum injection into the combined footprint
#                     footprint_injection = 0.0
#                     for linear_idx in combined_footprint
#                         cart_idx = CartesianIndices((nx, ny))[linear_idx]
#                         i_coord, j_coord = cart_idx.I
#                         footprint_injection += rain_rate_swim[i_coord, j_coord]
#                     end

#                     total_leaked += footprint_injection * effective_duration
#                 end
#             end
#         end

#         # Store the leaked volume for this trap
#         leaked_volumes_per_trap[trap_id] = total_leaked
#     end

#     return leaked_volumes_per_trap
# end

"""
Get total leaked volume (sum across all traps).

Parameters:
- leaked_volumes_per_trap: Dictionary from compute_per_trap_leaked_volumes()

Returns:
- Total leaked volume in swim units
"""
function get_total_leaked_volume(leaked_volumes_per_trap::Dict{Int, Float64})
    return sum(values(leaked_volumes_per_trap))
end

"""
Helper function to compute the height of CO2 in a specific trap at the current state.
Returns the maximum height within the trap's footprint.
"""
function _compute_trap_height(trap_id::Int, cur_amounts::Vector, z_vol_tables, tstruct)
    water_volume = cur_amounts[trap_id].amount

    # If trap is empty, height is zero
    if water_volume <= 0.0
        return 0.0
    end

    # Get water elevation from volume
    water_elevation = _volume_to_elevation(water_volume, z_vol_tables[trap_id])

    min_base_elevation = minimum(tstruct.topography[tstruct.footprints[trap_id]])
    distance_from_spillpoint_to_trap_bottom = 0.0
    children = subtrapsof(tstruct, trap_id)
    if !isempty(children)
        distance_from_spillpoint_to_trap_bottom = max.(distance_from_spillpoint_to_trap_bottom, tstruct.spillpoints[children[1]].elevation - min_base_elevation)
    end

    water_elevation += distance_from_spillpoint_to_trap_bottom

    # Height is elevation difference
    return max(0.0, water_elevation - min_base_elevation)
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

"""
Trigger leakage for a trap and propagate to all parent traps.
Modifies leakage_state and rain_rate in place.

The leaked volume is only recorded for the trap that triggered the leakage (trap_id).
Parent traps are marked as leaked to stop injection, but their volumes are not added
to the total leaked amount (to avoid double counting).

Uses source_tracker to determine which injection sources contributed CO2 to this trap,
and turns off injection only in those source regions.

Computes the total volume that will leak based on:
- Only filled parent traps connected to the leaking trap contribute
- Only (1 - sand_residual_co2_saturation) fraction of the volume leaks
- The remaining fraction stays as residually trapped CO2
"""
function _trigger_leakage!(trap_id::Int, leakage_state::LeakageState,
                          cur_amounts::Vector, rain_rate::Matrix{Float64},
                          z_vol_tables, tstruct, cur_time::Float64,
                          source_tracker::SourceTracker, filled_traps::Vector{Bool},
                          reservoir_properties::ReservoirProperties, verbose::Bool)
    # Get all traps that need to be marked as leaked (trap + all parents)
    all_affected_traps = _identify_all_parent_traps(trap_id, tstruct)

    # Record first leakage time if this is the first leakage
    if isnothing(leakage_state.first_leakage_time)
        leakage_state.first_leakage_time = cur_time
    end

    # Find which injection sources contributed to this leakage
    contributing_sources = find_regions_to_stop(source_tracker, trap_id)

    # First, record the leakage event for the trap that triggered it
    if !leakage_state.leaked_traps[trap_id]
        volume_in_trap = cur_amounts[trap_id].amount
        height_at_detection = _compute_trap_height(trap_id, cur_amounts, z_vol_tables, tstruct)

        # Find connected filled traps (these will contribute to leakage)
        connected_filled_traps = _find_connected_filled_parents(trap_id, filled_traps, tstruct)

        # Calculate total volume in connected filled traps
        total_connected_volume = 0.0
        for connected_trap in connected_filled_traps
            total_connected_volume += cur_amounts[connected_trap].amount
        end

        # Calculate leakable vs residual volumes
        # Only (1 - sand_residual_co2_saturation) fraction leaks through
        leakage_fraction = 1.0 - reservoir_properties.sand_residual_co2_saturation
        total_leakable_volume = total_connected_volume * leakage_fraction
        residual_trapped_volume = total_connected_volume * reservoir_properties.sand_residual_co2_saturation

        leakage_state.leaked_traps[trap_id] = true

        push!(leakage_state.leakage_events,
              LeakageEvent(cur_time, trap_id, height_at_detection, volume_in_trap,
                          copy(contributing_sources), connected_filled_traps,
                          total_leakable_volume, residual_trapped_volume))

        if verbose
            println("  Connected filled traps: $connected_filled_traps")
            println("  Total connected volume: $total_connected_volume")
            println("  Leakable volume: $total_leakable_volume ($(leakage_fraction*100)%)")
            println("  Residual trapped: $residual_trapped_volume ($(reservoir_properties.sand_residual_co2_saturation*100)%)")
        end
    end

    # Now mark all affected traps as leaked
    for affected_trap in all_affected_traps
        # Skip if already leaked
        if leakage_state.leaked_traps[affected_trap]
            continue
        end

        # Mark trap as leaked (but don't count volume again)
        leakage_state.leaked_traps[affected_trap] = true
    end

    if verbose
        println("  Turning off injection in source regions: $contributing_sources")
    end

    # Zero out injection in all regions that contributed to this leakage
    for source_region in contributing_sources
        region_mask = tstruct.regions .== source_region
        rain_rate[region_mask] .= 0.0
    end
end


function find_leakage_location(trap_id::Int, tstruct::TrapStructure)
    # Get the exact leakage location
    # This is the location in the trap footprint with the smallest topography elevation
    footprint = tstruct.footprints[trap_id]                # Vector of linear indices
    topography = tstruct.topography[footprint]             # Elevations at those indices
    min_idx = argmin(topography)                           # Index of minimum elevation within the footprint
    leakage_linear_idx = footprint[min_idx]                # Linear index in the domain
    leakage_cartesian_idx = CartesianIndices(size(tstruct.topography))[leakage_linear_idx]  # (i, j) coordinates
    return leakage_cartesian_idx
end


"""
Check all traps for leakage condition and trigger leakage if needed.
Returns true if any leakage was detected and triggered.
"""
function _check_and_trigger_leakage!(leakage_state::LeakageState, cur_amounts::Vector,
                                     rain_rate::Matrix{Float64}, z_vol_tables, tstruct,
                                     cur_time::Float64, source_tracker::SourceTracker,
                                     filled_traps::Vector{Bool}, reservoir_properties::ReservoirProperties,
                                     verbose::Bool)
    leakage_occurred = false
    num_traps = length(cur_amounts)

    for trap_id in 1:num_traps
        # Skip if already leaked
        if leakage_state.leaked_traps[trap_id]
            continue
        end

        # Compute current height
        height = _compute_trap_height(trap_id, cur_amounts, z_vol_tables, tstruct)

        # Check if exceeds threshold
        if height > leakage_state.leakage_threshold
            if verbose
                println("Leakage detected in trap $trap_id at time $cur_time (height: $height m)")
            end

            _trigger_leakage!(trap_id, leakage_state, cur_amounts, rain_rate,
                            z_vol_tables, tstruct, cur_time, source_tracker,
                            filled_traps, reservoir_properties, verbose)
            leakage_occurred = true
        end
    end

    return leakage_occurred
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
                                                       filled_traps, reservoir_properties, verbose)

        # If leakage occurred, we need to recompute flow rates with modified injection
        if leakage_occurred
            # Recompute flow with updated rain_rate (zeroed out in leaked regions)
            SurfaceWaterIntegratedModeling.setsavepoint!(rateinfo)
            rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, rain_rate,
                                                                   infiltration, tstruct, verbose)

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