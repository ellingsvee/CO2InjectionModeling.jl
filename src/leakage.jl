using SurfaceWaterIntegratedModeling
import Interpolations
using Graphs

export spillpoint_heights_above_trap_bottom
export get_loc_and_time_of_leakage

export get_leakage_root_trap
export get_leakage_trap


# function get_leakage_root_trap(
#     tstruct::TrapStructure{<:Real},
#     injection_location::Int,
#     leakage_height::Float64;
#     spillpoint_heights::Vector{Float64}
# )::Union{Int,Nothing}
#     N = numtraps(tstruct)

#     # Identify traps containing the injection location
#     # traps_in_footprint = findall(trap -> injection_location in tstruct.footprints[trap], 1:N)
#     traps_in_footprint = findall(trap -> injection_location in tstruct.regions[trap], 1:N)
#     if isempty(traps_in_footprint)
#         println("Injection location was ", injection_location)
#         println("Injection location not found in any trap footprint.")
#         return nothing
#     end

#     # Find traps where the spillpoint height is greater than leakage_height
#     traps_higher_than_leakage = findall(x -> x >= leakage_height, spillpoint_heights)
#     if isempty(traps_higher_than_leakage)
#         println("No traps have a spillpoint height >= $leakage_height.")
#         return nothing
#     end

#     candidate_traps = intersect(traps_higher_than_leakage, traps_in_footprint)

#     if isempty(candidate_traps)
#         println("No trap has a reachable leakage height for z_target.")
#         return nothing
#     end

#     # Select the trap with the smallest spillpoint height above trap bottom
#     root_trap = candidate_traps[argmin(spillpoint_heights[candidate_traps])]
#     return root_trap
# end

function all_subtraps(tstruct::TrapStructure{<:Real}, trap_ix::Int)
    subs = subtrapsof(tstruct, trap_ix)
    all_subs = copy(subs)
    for sub in subs
        append!(all_subs, all_subtraps(tstruct, sub))
    end
    return all_subs
end

# function get_leakage_trap(
#     tstruct::TrapStructure{<:Real},
#     leakage_height::Float64,
#     root_trap::Int;
#     zvtable::Vector{Tuple{Vector{Float64},Vector{Float64}}},
#     spillpoint_heights::Vector{Float64}
# )::Union{Int,Nothing}
#     # Identify traps contained within the footprint of the root trap
#     children = all_subtraps(tstruct, root_trap)
#     if !isempty(children)
#         traps_in_footprint = children
#         push!(traps_in_footprint, root_trap) # Include the root trap itself
#     else
#         traps_in_footprint = [root_trap]
#     end

#     # Find traps_in_footprint where the spillpoint height is greater than leakage_height
#     traps_in_footprint_higher_than_leakage = [trap for trap in traps_in_footprint if spillpoint_heights[trap] >= leakage_height]
#     if isempty(traps_in_footprint_higher_than_leakage)
#         println("No traps have a spillpoint height >= $leakage_height.")
#         return nothing
#     end

#     candidate_traps = traps_in_footprint_higher_than_leakage
#     sorted_candidates = sort(candidate_traps, by=x -> spillpoint_heights[x])

#     # Iterate through, and find the smallest trap that got zvals 
#     for trap in sorted_candidates
#         zvals = zvtable[trap][1]
#         zmin = minimum(zvals)
#         zmax = maximum(zvals)

#         topography_min = minimum(tstruct.topography[tstruct.footprints[trap]])
#         leakage_height_global = leakage_height + topography_min

#         if leakage_height_global >= zmin && leakage_height_global <= zmax
#             # println("Selected leakage trap: ", trap, "\n")
#             return trap
#         end
#     end

#     println("No trap has a reachable leakage height for z_target.")
#     return nothing
# end



function get_leakage_location(
    tstruct::TrapStructure{<:Real},
    leakage_root_trap::Int,
)
    # The leakage location will be the lowest point in the footprint of the leakage root trap
    footprint = tstruct.footprints[leakage_root_trap]
    leakage_location = footprint[argmin(tstruct.topography[footprint])]
    return leakage_location
end

function spillpoint_heights_above_trap_bottom(tstruct::TrapStructure{<:Real})
    N = numtraps(tstruct)

    heights = zeros(Float64, N)
    for i = 1:N
        spillpoint = tstruct.spillpoints[i]
        if !isnothing(spillpoint)
            trap_bottom = minimum(tstruct.topography[tstruct.footprints[i]])
            heights[i] = spillpoint.elevation - trap_bottom
        else
            heights[i] = 0.0
        end
    end

    return heights
end

"""
    time_to_reach_height(seq, tstruct, trap_ix, target_height; zvtable) -> Union{Float64, Nothing}

Calculate when a specific trap reaches a target height above its bottom.

Returns `nothing` if the target height is never reached or is outside the valid range.
"""
function time_to_reach_height(
    seq,
    tstruct::TrapStructure{<:Real},
    trap_ix::Int,
    target_height_above_trap_bottom::Float64;
    zvtable::Vector{Tuple{Vector{Float64},Vector{Float64}}}
)::Union{Float64,Nothing}
    zvals, vvals = zvtable[trap_ix]

    # Early return for invalid data
    if isempty(zvals) || isempty(vvals)
        return nothing
    end

    # Convert target height to global elevation
    trap_bottom = minimum(tstruct.topography[tstruct.footprints[trap_ix]])
    z_target = trap_bottom + target_height_above_trap_bottom

    # Check if target is within valid range
    z_min, z_max = extrema(zvals)
    if z_target < z_min || z_target > z_max
        return nothing
    end

    # Create interpolator for height-to-volume conversion
    z2v = if length(zvals) == 1
        _ -> vvals[1]  # Constant volume
    else
        Interpolations.linear_interpolation(
            zvals, vvals,
            extrapolation_bc=Interpolations.Flat()
        )
    end

    v_target = z2v(z_target)

    # Find the time when volume crosses the target
    for i in 2:length(seq)
        v_prev = amount_at(seq, i - 1)[trap_ix].amount
        v_curr = amount_at(seq, i)[trap_ix].amount

        # Check if target volume is crossed in this interval
        if (v_prev <= v_target <= v_curr) || (v_curr <= v_target <= v_prev)
            t_prev = seq[i-1].timestamp
            t_curr = seq[i].timestamp

            # Linear interpolation for crossing time
            if v_curr ≈ v_prev
                return t_prev  # Volume not changing, use start time
            end

            fraction = (v_target - v_prev) / (v_curr - v_prev)
            return t_prev + fraction * (t_curr - t_prev)
        end
    end

    return nothing  # Target not reached
end


"""
    fill_until_height(seq, tstruct, target_height; zvtable) -> Union{Tuple{Int, Float64}, Nothing}

Find which trap reaches the target height first and when.

Returns `(trap_index, time)` for the earliest trap to reach the height,
or `nothing` if no trap reaches it.
"""
function fill_until_height(
    seq,
    tstruct::TrapStructure{<:Real},
    target_height_above_trap_bottom::Float64;
    zvtable::Vector{Tuple{Vector{Float64},Vector{Float64}}}
)::Union{Tuple{Int,Float64},Nothing}
    N = numtraps(tstruct)

    # Compute fill times for all traps
    fill_times = [
        time_to_reach_height(seq, tstruct, trap_ix, target_height_above_trap_bottom; zvtable)
        for trap_ix in 1:N
    ]

    # Find the earliest valid fill time
    valid_indices = findall(!isnothing, fill_times)

    if isempty(valid_indices)
        return nothing
    end

    # Get trap with minimum fill time
    earliest_idx = argmin(fill_times[valid_indices])
    earliest_trap = valid_indices[earliest_idx]
    earliest_time = fill_times[earliest_trap]

    return (earliest_trap, earliest_time)
end

function get_loc_and_time_of_leakage(seq, tstruct, injection_location, leakage_height)
    zvtable = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # Use the simpler approach: find first trap that reaches target height
    result = fill_until_height(seq, tstruct, leakage_height; zvtable=zvtable)

    if isnothing(result)
        println("No traps reached the target height - CO2 has exited the domain.")
        return (nothing, nothing)
    end

    earliest_trap, earliest_time = result
    leakage_location = get_leakage_location(tstruct, earliest_trap)

    return (leakage_location, earliest_time)
end
