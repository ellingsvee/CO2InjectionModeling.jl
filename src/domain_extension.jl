# Utility function to extend domain boundaries with a dome-like structure
export extend_domain_with_dome

"""
    extend_domain_with_dome(topographies, boundary_width, target_height; transition_type=:exponential)

Extend the domain boundaries to create a dome-like structure that prevents CO2 from flowing out
of the domain while maintaining natural-looking geometry.

# Arguments
- `topographies::Array{Float64,3}`: Original topography data (n_layers × ny × nx)
- `boundary_width::Int`: Number of cells to add on each side as boundary
- `target_height::Union{Float64,Nothing}=nothing`: Target height for the dome boundary.
    If `nothing`, uses `maximum(topographies[i, :, :]) + 10` for each layer.
- `transition_type::Symbol=:exponential`: Type of transition (:linear, :exponential, :sigmoid)
- `steepness::Float64=2.0`: Controls how quickly the boundary rises (higher = steeper)

# Returns
- `extended_topographies::Array{Float64,3}`: Extended topography with dome boundaries
- `inner_indices::Tuple`: Indices of the original (inner) domain as `(i_range, j_range)`

# Example
```julia
# Extend domain with 50-cell boundary that rises exponentially
extended_topo, (i_range, j_range) = extend_domain_with_dome(
    caprock_topography, 50; 
    transition_type=:exponential, 
    steepness=2.0
)

# Use extended_topo for simulation
tstructs = analyze_layers(extended_topo; lengths=(length_x, length_y))

# For visualization, use only the inner domain:
# inner_topo = extended_topo[:, i_range, j_range]
```
"""
function extend_domain_with_dome(
    topographies::Array{Float64,3},
    boundary_width::Int;
    target_height::Union{Float64,Nothing}=nothing,
    transition_type::Symbol=:exponential,
    steepness::Float64=2.0
)
    n_layers, ny, nx = size(topographies)
    
    # New dimensions
    new_ny = ny + 2 * boundary_width
    new_nx = nx + 2 * boundary_width
    
    # Initialize extended array
    extended_topographies = zeros(Float64, n_layers, new_ny, new_nx)
    
    # Process each layer
    for layer_idx in 1:n_layers
        original_topo = topographies[layer_idx, :, :]
        
        # Determine target height for this layer
        if isnothing(target_height)
            layer_target_height = maximum(original_topo) + 10.0
        else
            layer_target_height = target_height
        end
        
        # Copy original data to center
        extended_topographies[layer_idx, (boundary_width+1):(boundary_width+ny), (boundary_width+1):(boundary_width+nx)] .= original_topo
        
        # Create distance field from edge
        for i in 1:new_ny
            for j in 1:new_nx
                # Skip if in original domain
                if i > boundary_width && i <= boundary_width + ny && 
                   j > boundary_width && j <= boundary_width + nx
                    continue
                end
                
                # Calculate distance from boundary edge (how far into the boundary region)
                dist_from_edge = 0
                
                # Distance from top/bottom edges
                if i <= boundary_width
                    dist_from_edge = max(dist_from_edge, boundary_width - i + 1)
                elseif i > boundary_width + ny
                    dist_from_edge = max(dist_from_edge, i - (boundary_width + ny))
                end
                
                # Distance from left/right edges
                if j <= boundary_width
                    dist_from_edge = max(dist_from_edge, boundary_width - j + 1)
                elseif j > boundary_width + nx
                    dist_from_edge = max(dist_from_edge, j - (boundary_width + nx))
                end
                
                # Get nearest edge value from original domain
                nearest_i = clamp(i - boundary_width, 1, ny)
                nearest_j = clamp(j - boundary_width, 1, nx)
                edge_value = original_topo[nearest_i, nearest_j]
                
                # Calculate transition factor (0 at original edge, 1 at boundary edge)
                t = dist_from_edge / boundary_width
                
                # Apply transition function
                if transition_type == :linear
                    blend_factor = t
                elseif transition_type == :exponential
                    blend_factor = (exp(steepness * t) - 1) / (exp(steepness) - 1)
                elseif transition_type == :sigmoid
                    # Sigmoid transition centered at t=0.5
                    x = (t - 0.5) * steepness * 2
                    blend_factor = 1 / (1 + exp(-x))
                else
                    error("Unknown transition_type: $transition_type. Use :linear, :exponential, or :sigmoid")
                end
                
                # Interpolate between edge value and target height
                extended_topographies[layer_idx, i, j] = edge_value + blend_factor * (layer_target_height - edge_value)
            end
        end
    end
    
    # Return extended topographies and indices for original domain
    inner_indices = ((boundary_width+1):(boundary_width+ny), (boundary_width+1):(boundary_width+nx))
    
    return extended_topographies, inner_indices
end
