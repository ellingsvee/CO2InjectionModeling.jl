using SurfaceWaterIntegratedModeling
export SourceTracker, update_sources_on_trap_fill!, find_contributing_sources, find_regions_to_stop, update_injection_sources!

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