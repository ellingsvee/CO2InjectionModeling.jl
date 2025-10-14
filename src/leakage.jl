using SurfaceWaterIntegratedModeling
import Interpolations
using Graphs

export spillpoint_heights_above_trap_bottom
export get_loc_and_time_of_leakage

export get_leakage_root_trap
export get_leakage_trap


function get_leakage_root_trap(
    tstruct::TrapStructure{<:Real},
    injection_location::Int,
    leakage_height::Float64;
    spillpoint_heights::Vector{Float64}
)::Union{Int,Nothing}
    N = numtraps(tstruct)

    # Identify traps containing the injection location
    traps_in_footprint = findall(trap -> injection_location in tstruct.footprints[trap], 1:N)
    if isempty(traps_in_footprint)
        println("Injection location not found in any trap footprint.")
        return nothing
    end

    # Find traps where the spillpoint height is greater than leakage_height
    traps_higher_than_leakage = findall(x -> x >= leakage_height, spillpoint_heights)
    if isempty(traps_higher_than_leakage)
        println("No traps have a spillpoint height >= $leakage_height.")
        return nothing
    end

    candidate_traps = intersect(traps_higher_than_leakage, traps_in_footprint)

    if isempty(candidate_traps)
        println("No trap has a reachable leakage height for z_target.")
        return nothing
    end

    # Select the trap with the smallest spillpoint height above trap bottom
    root_trap = candidate_traps[argmin(spillpoint_heights[candidate_traps])]
    return root_trap
end

function all_subtraps(tstruct::TrapStructure{<:Real}, trap_ix::Int)
    subs = subtrapsof(tstruct, trap_ix)
    all_subs = copy(subs)
    for sub in subs
        append!(all_subs, all_subtraps(tstruct, sub))
    end
    return all_subs
end

function get_leakage_trap(
    tstruct::TrapStructure{<:Real},
    leakage_height::Float64,
    root_trap::Int;
    zvtable::Vector{Tuple{Vector{Float64},Vector{Float64}}},
    spillpoint_heights::Vector{Float64}
)::Union{Int,Nothing}
    # Identify traps contained within the footprint of the root trap
    children = all_subtraps(tstruct, root_trap)
    if !isempty(children)
        traps_in_footprint = children
        push!(traps_in_footprint, root_trap) # Include the root trap itself
    else
        traps_in_footprint = [root_trap]
    end

    # Find traps_in_footprint where the spillpoint height is greater than leakage_height
    traps_in_footprint_higher_than_leakage = [trap for trap in traps_in_footprint if spillpoint_heights[trap] >= leakage_height]
    if isempty(traps_in_footprint_higher_than_leakage)
        println("No traps have a spillpoint height >= $leakage_height.")
        return nothing
    end

    # println("candidate_traps = ", traps_in_footprint_higher_than_leakage, "\n")
    candidate_traps = traps_in_footprint_higher_than_leakage
    sorted_candidates = sort(candidate_traps, by=x -> spillpoint_heights[x])

    # Iterate through, and find the smallest trap that got zvals 
    # println("length(sorted_candidates) = ", length(sorted_candidates), "\n")
    for trap in sorted_candidates
        zvals = zvtable[trap][1]
        zmin = minimum(zvals)
        zmax = maximum(zvals)

        topography_min = minimum(tstruct.topography[tstruct.footprints[trap]])
        leakage_height_global = leakage_height + topography_min

        # println("Trap: ", trap, ", zmin: ", zmin, ", zmax: ", zmax, ", topography_min: ", topography_min, ", leakage_height_global: ", leakage_height_global, "\n")

        if leakage_height_global >= zmin && leakage_height_global <= zmax
            # println("Selected leakage trap: ", trap, "\n")
            return trap
        end
    end

    println("No trap has a reachable leakage height for z_target.")
    return nothing
end



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

function time_to_reach_height(seq, tstruct, trap_ix, target_height_above_trap_bottom; zvtable::Vector{Tuple{Vector{Float64},Vector{Float64}}})
    # zvtable = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    zvals, vvals = zvtable[trap_ix]

    # Compute trap bottom consistently
    trap_bottom_global = minimum(tstruct.topography[tstruct.footprints[trap_ix]])
    z_target = trap_bottom_global + target_height_above_trap_bottom

    # Check if z_target is within zvals range
    if z_target < minimum(zvals) || z_target > maximum(zvals)
        println("z_target ($z_target) is outside zvals range [$(minimum(zvals)), $(maximum(zvals))]")
        return nothing
    end

    z2v = length(zvals) == 1 ?
          (z -> 0.0) :
          Interpolations.linear_interpolation(zvals, vvals, extrapolation_bc=Interpolations.Flat())

    v_target = z2v(z_target)

    for i in 2:length(seq)
        v_prev = amount_at(seq, i - 1)[trap_ix].amount
        v_next = amount_at(seq, i)[trap_ix].amount
        t_prev = seq[i-1].timestamp
        t_next = seq[i].timestamp
        if (v_prev < v_target && v_next >= v_target) || (v_prev > v_target && v_next <= v_target)
            frac = (v_target - v_prev) / (v_next - v_prev)
            t_cross = t_prev + frac * (t_next - t_prev)
            return t_cross
        end
    end
    println("Target height not reached in the sequence.")
    return nothing
end

function get_loc_and_time_of_leakage(seq, tstruct, injection_location, leakage_height)
    zvtable = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)
    spillpoint_heights = spillpoint_heights_above_trap_bottom(tstruct)

    # Find the trap that is leaking
    root_trap = get_leakage_root_trap(tstruct, injection_location, leakage_height, spillpoint_heights=spillpoint_heights)

    if isnothing(root_trap)
        println("No leaking trap found.")
        return nothing
    else
        println("Root trap: ", root_trap)
    end

    leakage_trap = get_leakage_trap(tstruct, leakage_height, root_trap, zvtable=zvtable, spillpoint_heights=spillpoint_heights)

    if isnothing(leakage_trap)
        println("No leakage trap found.")
        return nothing
    else
        println("Leakage trap: ", leakage_trap)
    end

    # Get the time at which the leakage occurs
    time_at_leakage = time_to_reach_height(seq, tstruct, leakage_trap, leakage_height; zvtable=zvtable)
    leakage_location = get_leakage_location(tstruct, leakage_trap)
    return (leakage_location, time_at_leakage)
end