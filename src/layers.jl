using SurfaceWaterIntegratedModeling

export analyze_layers
export run_injection

"""
    analyze_layers(layers; lengths=nothing)

Analyze trap structures for multiple layers.

# Arguments
- `layers::Vector{Matrix{Float64}}` or `Array{Float64,3}`: Topography grids for each layer
- `lengths::Union{Tuple{<:Real, <:Real}, Nothing}=nothing`: Physical dimensions (length_x, length_y) of the grid in meters

# Returns
- `Vector{TrapStructure{Float64}}`: Trap structures for each layer

# Example
```julia
# Without physical dimensions (grid units assumed)
tstructs = analyze_layers(layers)

# With physical dimensions from Sleipner data
# If grid is 461×701 cells covering (xmin=435000, xmax=470000) and (ymin=6475000, ymax=6498000)
length_x = 470000 - 435000  # 35000 meters
length_y = 6498000 - 6475000  # 23000 meters
tstructs = analyze_layers(layers; lengths=(length_x, length_y))
```
"""
# function analyze_layers(
#     layers::Vector{Matrix{Float64}};
# )::Vector{TrapStructure{Float64}}
#     tstructs = TrapStructure{Float64}[]
#     for layer in layers
#         push!(tstructs, spillanalysis(layer))
#     end
#     return tstructs
# end

# function analyze_layers(
#     layers::Array{Float64,3};
# )::Vector{TrapStructure{Float64}}
#     tstructs = TrapStructure{Float64}[]

#     for i in 1:size(layers, 1)
#         push!(tstructs, spillanalysis(layers[i, :, :]))
#     end
#     return tstructs
# end
function analyze_layers(
    layers::Vector{Matrix{Float64}};
    lengths::Union{Tuple{<:Real,<:Real},Nothing}=nothing
)::Vector{TrapStructure{Float64}}
    tstructs = TrapStructure{Float64}[]
    for layer in layers
        push!(tstructs, spillanalysis(layer; lengths=lengths))
    end
    return tstructs
end

function analyze_layers(
    layers::Array{Float64,3};
    lengths::Union{Tuple{<:Real,<:Real},Nothing}=nothing
)::Vector{TrapStructure{Float64}}
    tstructs = TrapStructure{Float64}[]
    for i in 1:size(layers, 1)
        push!(tstructs, spillanalysis(layers[i, :, :]; lengths=lengths))
        # push!(tstructs, spillanalysis(layers[i, :, :]))
    end
    return tstructs
end

"""
Function that runs a full injection simulation with leakage between the layers.
Takes in the tstructs, injection location in the top layer, the injection rate, leakage heights for the different layers.
"""
function run_injection(
    tstructs::Vector{TrapStructure{Float64}},
    injection_location::Int,
    injection_rate::Float64,
    leakage_heights::Vector{Float64}
)
    @assert length(tstructs) == length(leakage_heights) + 1 "Number of leakage heights must be one less than number of layers"


    topography_top_layer = tstructs[1].topography
    injection_rate_matrix = fill(0.0, size(topography_top_layer))
    injection_rate_matrix[injection_location] = injection_rate
    injection = [WeatherEvent(0.0, injection_rate_matrix)] # Register the injection as a WeatherEvent

    seqs = Vector{Vector{SpillEvent}}(undef, length(tstructs))

    leakage_locations = Vector{Union{Nothing,Int}}(undef, length(tstructs) - 1)
    times_at_leakage = Vector{Union{Nothing,Float64}}(undef, length(tstructs) - 1)


    # Initialize leakage_location and time_at_leakage for the first layer
    leakage_location = injection_location
    time_at_leakage = 0.0

    for (i, tstruct) in enumerate(tstructs)
        println("Simulating layer $i...")

        # Register the injection/leakage as a WeatherEvent
        topography = tstructs[i].topography
        injection_rate_matrix = fill(0.0, size(topography))
        injection_rate_matrix[leakage_location] = injection_rate
        # injection = [WeatherEvent(time_at_leakage, injection_rate_matrix)]

        # Define the weather sequence with the current injection/leakage
        if i > 1
            # If not the first layer, ensure no injection until the leakage time
            injection = [
                WeatherEvent(0.0, fill(0.0, size(topography))), # No injection until the leakage time
                WeatherEvent(time_at_leakage, injection_rate_matrix)
            ]
        else
            # For the first layer, start injection at time 0
            injection = [WeatherEvent(0.0, injection_rate_matrix)]
        end

        # Run the fill sequence
        seqs[i] = fill_sequence(tstruct, injection)

        # No leakage for the last layer
        if i == length(tstructs)
            # Use the last timestamp in the sequence for consistency
            # times_at_leakage[i] = seqs[i][end].timestamp
            # leakage_locations[i] = nothing
            break
        end

        # Analyzing to get the time and location of leakage
        leakage_height = leakage_heights[i]
        leakage_result = get_loc_and_time_of_leakage(seqs[i], tstruct, leakage_location, leakage_height)

        # Check if leakage occurred
        if isnothing(leakage_result) || isnothing(leakage_result[1]) || isnothing(leakage_result[2])
            println("No leakage occurred in layer $i - CO2 has likely exited the domain.")
            leakage_locations[i] = nothing
            times_at_leakage[i] = nothing

            # For remaining layers, run fill_sequence with zero injection
            # This maintains consistent sequence lengths and structure
            for j in (i+1):length(tstructs)
                println("Layer $j: Running with zero injection (CO2 already exited domain)")
                topography_j = tstructs[j].topography
                zero_injection = [WeatherEvent(0.0, fill(0.0, size(topography_j)))]
                seqs[j] = fill_sequence(tstructs[j], zero_injection)

                # Mark as no leakage for remaining layers
                if j < length(tstructs)
                    leakage_locations[j] = nothing
                    times_at_leakage[j] = nothing
                end
            end
            break
        else
            # Leakage occurred successfully
            (leakage_location, time_at_leakage) = leakage_result
            leakage_locations[i] = leakage_location
            times_at_leakage[i] = time_at_leakage
            println("Leakage occurred in layer $i at time $time_at_leakage at location $leakage_location.")
        end
    end

    # start_time = 0.0
    # end_time = seqs[end][end].timestamp

    # Generate times with start_time, times_at_leakage, and end_time
    # times = [start_time + 0.1; [t for t in times_at_leakage if t !== nothing]; end_time] # Avoid exact start_time to prevent interpolation issues

    return seqs, times_at_leakage, leakage_locations
end