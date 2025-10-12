using SurfaceWaterIntegratedModeling
import Interpolations

export spillpoint_heights_above_trap_bottom
export get_time_of_leakage

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

function time_to_reach_height(seq, tstruct, trap_ix, target_height_above_trap_bottom)
    zvtable = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)
    zvals, vvals = zvtable[trap_ix]

    # Use the global bottom (lowest point in the footprint)
    trap_bottom_global = minimum(tstruct.topography[tstruct.footprints[trap_ix]])

    # Target height above the global bottom
    z_target = trap_bottom_global + target_height_above_trap_bottom

    z2v = length(zvals) == 1 ?
          (z -> 0.0) :
          Interpolations.linear_interpolation(zvals, vvals, extrapolation_bc=Interpolations.Line())

    v_target = z2v(z_target)
    for i in 2:length(seq)
        v_prev = amount_at(seq, i - 1)[trap_ix].amount
        v_next = amount_at(seq, i)[trap_ix].amount
        t_prev = seq[i-1].timestamp
        t_next = seq[i].timestamp
        if (v_prev < v_target && v_next >= v_target) || (v_prev > v_target && v_next <= v_target)
            frac = (v_target - v_prev) / (v_next - v_prev)
            t_cross = t_prev + frac * (t_next - t_prev)
            # return t_cross, v_target
            return t_cross
        end
    end
    println("Target height not reached in the sequence.")
    return nothing
end

function get_time_of_leakage(seq, tstruct, injection_location, leakage_height)
    # Find the trap that is leaking
    leaking_trap = get_leakage_root_trap(tstruct, injection_location, leakage_height)
    if isnothing(leaking_trap)
        println("No leaking trap found.")
        return nothing
    end

    # Get the time at which the leakage occurs
    time_at_leakage = time_to_reach_height(seq, tstruct, leaking_trap, leakage_height)
    return time_at_leakage
end