export two_dim_index_to_linear_index
export linear_index_to_two_dim_index
export compute_co2_volumes_over_time
export compute_ensemble_volume_statistics
export compute_injection_metrics


function two_dim_index_to_linear_index(grid, idxs)
    return LinearIndices(grid)[idxs...]
end

function linear_index_to_two_dim_index(grid, idx)
    return Tuple(CartesianIndices(grid)[idx])
end

"""
    compute_co2_volumes_over_time(tstructs, all_texs, common_times, cell_volume)

Compute the volume of trapped CO2 in each layer over time for all simulations.

# Arguments
- `tstructs::Vector{TrapStructure}`: Trap structures for all layers
- `all_texs::Vector{Vector{Vector}}`: Textures for each simulation, layer, and time
- `common_times::Vector{Float64}`: Time points
- `cell_volume::Float64`: Volume of a single grid cell (length_x * length_y / (nx * ny))

# Returns
- `all_volumes::Vector{Vector{Vector{Float64}}}`: Volumes for each simulation, layer, and time
  - `all_volumes[sim][layer][time_idx]` = volume of CO2 in that layer at that time
"""
function compute_co2_volumes_over_time(
    tstructs::Vector,
    all_texs::Vector,
    common_times::Vector{Float64},
    cell_volume::Float64
)
    n_simulations = length(all_texs)
    n_layers = length(tstructs)
    n_times = length(common_times)
    
    all_volumes = Vector{Vector{Vector{Float64}}}(undef, n_simulations)
    
    for sim_idx in 1:n_simulations
        layer_volumes = Vector{Vector{Float64}}(undef, n_layers)
        
        for layer_idx in 1:n_layers
            time_volumes = Vector{Float64}(undef, n_times)
            
            for time_idx in 1:n_times
                # Get texture at this time
                tex = all_texs[sim_idx][layer_idx][time_idx]
                
                # Count cells with CO2 (tex > 1 indicates CO2 presence)
                n_filled_cells = count(tex .> 1)
                
                # Compute volume (number of cells × cell volume)
                time_volumes[time_idx] = n_filled_cells * cell_volume
            end
            
            layer_volumes[layer_idx] = time_volumes
        end
        
        all_volumes[sim_idx] = layer_volumes
    end
    
    return all_volumes
end

"""
    compute_ensemble_volume_statistics(all_volumes)

Compute mean and standard deviation of CO2 volumes across ensemble simulations.

# Arguments
- `all_volumes::Vector{Vector{Vector{Float64}}}`: Volumes from `compute_co2_volumes_over_time`

# Returns
- `mean_volumes::Vector{Vector{Float64}}`: Mean volume for each layer and time
- `std_volumes::Vector{Vector{Float64}}`: Standard deviation for each layer and time
- `max_volumes::Vector{Float64}`: Maximum possible volume for each layer (from trap capacity)
"""
function compute_ensemble_volume_statistics(all_volumes::Vector{Vector{Vector{Float64}}})
    n_simulations = length(all_volumes)
    n_layers = length(all_volumes[1])
    n_times = length(all_volumes[1][1])
    
    mean_volumes = Vector{Vector{Float64}}(undef, n_layers)
    std_volumes = Vector{Vector{Float64}}(undef, n_layers)
    
    for layer_idx in 1:n_layers
        mean_vol = Vector{Float64}(undef, n_times)
        std_vol = Vector{Float64}(undef, n_times)
        
        for time_idx in 1:n_times
            # Collect volumes from all simulations for this layer and time
            volumes_at_time = [all_volumes[sim][layer_idx][time_idx] for sim in 1:n_simulations]
            
            mean_vol[time_idx] = mean(volumes_at_time)
            std_vol[time_idx] = std(volumes_at_time)
        end
        
        mean_volumes[layer_idx] = mean_vol
        std_volumes[layer_idx] = std_vol
    end
    
    # Compute max volumes (maximum observed volume across all simulations and times)
    max_volumes = Vector{Float64}(undef, n_layers)
    for layer_idx in 1:n_layers
        max_vol = 0.0
        for sim_idx in 1:n_simulations
            max_vol = max(max_vol, maximum(all_volumes[sim_idx][layer_idx]))
        end
        max_volumes[layer_idx] = max_vol
    end
    
    return mean_volumes, std_volumes, max_volumes
end

"""
    compute_injection_metrics(common_times, all_volumes, injection_rate, co2_density)

Compute injection-related metrics such as cumulative injected mass and storage efficiency.

# Arguments
- `common_times::Vector{Float64}`: Time points (years)
- `all_volumes::Vector{Vector{Vector{Float64}}}`: Volumes from `compute_co2_volumes_over_time`
- `injection_rate::Float64`: Injection rate (kg/s or specified units)
- `co2_density::Float64`: CO2 density (kg/m³) under storage conditions

# Returns
- `cumulative_injected::Vector{Float64}`: Cumulative mass injected over time (same units as rate×time)
- `mean_stored_fraction::Vector{Float64}`: Mean fraction of injected CO2 that is stored (0-1)
- `std_stored_fraction::Vector{Float64}`: Standard deviation of storage fraction
"""
function compute_injection_metrics(
    common_times::Vector{Float64},
    all_volumes::Vector{Vector{Vector{Float64}}},
    injection_rate::Float64,
    co2_density::Float64
)
    n_simulations = length(all_volumes)
    n_layers = length(all_volumes[1])
    n_times = length(common_times)
    
    # Convert years to seconds for rate calculation
    seconds_per_year = 365.25 * 24 * 3600
    
    # Cumulative injected mass (assuming constant injection rate)
    cumulative_injected = injection_rate .* (common_times .* seconds_per_year)
    
    # Compute stored mass for each simulation
    stored_fractions = Matrix{Float64}(undef, n_simulations, n_times)
    
    for sim_idx in 1:n_simulations
        for time_idx in 1:n_times
            # Total stored volume across all layers
            total_volume = sum([all_volumes[sim_idx][layer_idx][time_idx] for layer_idx in 1:n_layers])
            
            # Total stored mass
            stored_mass = total_volume * co2_density
            
            # Fraction of injected mass that is stored
            if cumulative_injected[time_idx] > 0
                stored_fractions[sim_idx, time_idx] = stored_mass / cumulative_injected[time_idx]
            else
                stored_fractions[sim_idx, time_idx] = 0.0
            end
        end
    end
    
    # Compute mean and std of storage fraction
    mean_stored_fraction = vec(mean(stored_fractions, dims=1))
    std_stored_fraction = vec(std(stored_fractions, dims=1))
    
    return cumulative_injected, mean_stored_fraction, std_stored_fraction
end
