import SurfaceWaterIntegratedModeling
import Interpolations

export compute_co2_height_matrix, create_co2_mask_3d_from_heights, remove_boundary_padding

"""
Compute a 2D matrix (or vector of matrices) of CO2 height above the layer at given time(s).

This function can work in multiple modes:
1. Single event index (idx::Int): Returns one height matrix from that event
2. Multiple event indices (idx::Vector{Int}): Returns vector of height matrices
3. Single timepoint (timepoint::Float64): Returns one interpolated height matrix
4. Multiple timepoints (timepoint::Vector{Float64}): Returns vector of interpolated height matrices

Parameters:
- seq: Vector of SpillEvents from the simulation
- tstruct: TrapStructure from the layer
- idx: (Optional) Index or vector of indices in seq to compute height matrices for.
       If not provided, timepoint must be specified.
- timepoint: (Optional) Single time value or vector of time values to compute height at.
            If provided, this overrides idx and uses interpolation.
- use_layer_base: (Optional, default=false) If true, heights are measured from the actual
                 layer topography. If false, heights are measured from the trap's effective
                 bottom (which for parent traps is above their subtraps, as defined in
                 _compute_z_vol_tables). Setting to true gives the total CO2 column height.

Returns:
- If single idx/timepoint: 2D Matrix{Float64} with CO2 heights
- If vector of idx/timepoints: Vector{Matrix{Float64}} with one matrix per time

Each height_matrix contains the height of CO2 above the reference level at each location.
Zero values indicate no CO2 is present.

Examples:
```julia
# Single event index - height from trap bottom
height = compute_co2_height_matrix(seq, tstruct, idx=5)

# Height from actual layer base (total column height)
height = compute_co2_height_matrix(seq, tstruct, idx=5, use_layer_base=true)

# Multiple event indices
heights = compute_co2_height_matrix(seq, tstruct, idx=[5, 10, 15])

# Single arbitrary timepoint (with interpolation)
height = compute_co2_height_matrix(seq, tstruct, timepoint=2.5)

# Multiple arbitrary timepoints (with interpolation)
heights = compute_co2_height_matrix(seq, tstruct, timepoint=[1.0, 2.5, 5.0])
```
"""
function compute_co2_height_matrix(
    seq::Vector{SurfaceWaterIntegratedModeling.SpillEvent},
    tstruct::TrapStructure;
    idx::Union{Int, Vector{Int}, Nothing} = nothing,
    timepoint::Union{Float64, Vector{Float64}, Nothing} = nothing,
    use_layer_base::Bool = false,
    tstates::Union{Nothing, Vector{Any}} = nothing,
)
    # Validate inputs
    if isnothing(idx) && isnothing(timepoint)
        error("Either idx or timepoint must be provided")
    end

    if !isnothing(idx) && !isnothing(timepoint)
        error("Cannot specify both idx and timepoint")
    end

    # Determine if we're computing for single or multiple times
    is_vector_input = false
    times_to_compute = nothing

    if !isnothing(timepoint)
        if timepoint isa Vector
            is_vector_input = true
            times_to_compute = timepoint
        else
            times_to_compute = [timepoint]
        end
    elseif !isnothing(idx)
        if idx isa Vector
            is_vector_input = true
            times_to_compute = nothing  # Will use indices directly
        else
            times_to_compute = nothing  # Will use single index
        end
    end

    # Compute z_vol_tables once (expensive operation)
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # Get topography dimensions
    nx, ny = size(tstruct.topography)

    # Helper function to compute height matrix for a single set of amounts
    function _compute_single_height_matrix(amounts::Vector)
        height_matrix = zeros(Float64, nx, ny)

        # Process traps in order from lowest-level (regions) to highest (parent traps)
        # This ensures child traps are processed before their parents
        for trap_id in 1:numtraps(tstruct)
            water_volume = amounts[trap_id].amount

            # Skip if trap is empty
            if water_volume <= 0.0
                continue
            end

            # Get the water level (elevation) for this trap given its volume
            # Note: The z_vol_tables and water_volume already account for the correct
            # trap geometry (parent traps have their bottom above children)
            water_elevation = _volume_to_elevation(water_volume, z_vol_tables[trap_id])

            # Get the footprint of this trap (cells belonging to it)
            footprint = tstruct.footprints[trap_id]

            # Determine the reference elevation for height calculation
            if use_layer_base
                # Use actual layer topography - gives total column height
                reference_elevation = tstruct.topography
            else
                # Use effective trap bottom (parent traps start above their subtraps)
                # This matches the logic in _compute_z_vol_tables
                reference_elevation = copy(tstruct.topography)
                children = subtrapsof(tstruct, trap_id)
                if !isempty(children)
                    # This is a parent trap - elevate bottom to child spillpoint
                    child_spillpoint_elevation = tstruct.spillpoints[children[1]].elevation
                    for linear_idx in footprint
                        cart_idx = CartesianIndices((nx, ny))[linear_idx]
                        i, j = cart_idx.I
                        reference_elevation[i, j] = max(reference_elevation[i, j],
                                                       child_spillpoint_elevation)
                    end
                end
            end

            # For each cell in the footprint, compute height above reference
            for linear_idx in footprint
                # Convert linear index to 2D coordinates
                cart_idx = CartesianIndices((nx, ny))[linear_idx]
                i, j = cart_idx.I

                # Reference elevation at this cell
                base_elevation = reference_elevation[i, j]

                # Height of CO2 above reference
                # Using max() handles overlapping traps correctly - we want the highest water level
                if water_elevation > base_elevation
                    height_matrix[i, j] = max(height_matrix[i, j],
                                             water_elevation - base_elevation)
                end
            end
        end

        return height_matrix
    end

    # Compute based on input type
    if !isnothing(timepoint)
        # Using timepoint(s) - need interpolation 
        if isnothing(tstates)
            # Precompute trap states at requested timepoints
            tstates = trap_states_at_timepoints(tstruct, seq, times_to_compute, verbose=false)
        end
        height_matrices = [
            _compute_single_height_matrix([
                SurfaceWaterIntegratedModeling.FilledAmount(amt, times_to_compute[i])
                for amt in tstates[i][2]
            ])
            for i in eachindex(times_to_compute, tstates)
        ]

        return is_vector_input ? height_matrices : height_matrices[1]

    else
        # Using idx/indices - direct lookup
        if idx isa Int
            amounts = amount_at(seq, idx)
            return _compute_single_height_matrix(amounts)
        else
            # Multiple indices
            height_matrices = [
                _compute_single_height_matrix(amount_at(seq, i))
                for i in idx
            ]
            return height_matrices
        end
    end
end

"""
Helper function to convert volume to elevation using the z-vol table.
Uses Interpolations.jl for robust interpolation, following the same pattern as swim.
"""
function _volume_to_elevation(volume::Float64, z_vol_table)
    # z_vol_table is a tuple: (z_values, volume_values)
    # We need to create an interpolator from volume -> elevation

    zvals = z_vol_table[1]  # Elevations
    vvals = z_vol_table[2]  # Volumes

    # Handle degenerate case (single point)
    if length(vvals) == 1
        return zvals[1]
    end

    # Create interpolator for volume-to-elevation conversion
    # This inverts the relationship: instead of z->v, we do v->z
    v2z = Interpolations.linear_interpolation(
        vvals, zvals,
        extrapolation_bc=Interpolations.Line()
    )

    return v2z(volume)
end

"""
Create a 3D mask of CO2-filled cells based on computed height matrix.

This function converts a 2D height matrix (from compute_co2_height_matrix) into a 3D
boolean mask indicating which cells in the domain are filled with CO2. This provides
a more accurate representation than assuming traps are completely filled.

The mask always uses the layer's base topography as the reference. The CO2 fills from
the base topography upward by the height specified in height_matrix at each (i,j) location.

Parameters:
- height_matrix: 2D Matrix{Float64} from compute_co2_height_matrix, containing CO2 heights
                above the layer base at each (x,y) location. Note: If computed with
                use_layer_base=false, the heights represent CO2 in trap volumes only.
                For total column height from base, use use_layer_base=true when computing.
                The height_matrix may be padded if the layer uses closed boundary conditions.
- layer: Layer struct containing the trap_structure with topography information
- domain: Domain3D struct defining the 3D grid (original unpadded dimensions)

Returns:
- mask_3d: 3D boolean array (nx, ny, nz) where true indicates cells filled with CO2

The mask is computed by:
1. For each (i,j) horizontal location with non-zero height
2. Find all vertical cells (k) where the cell depth is between the base and (base + height)
3. Mark those cells as filled

Examples:
```julia
# Compute height matrix at a specific time (with total column height)
height_matrix = compute_co2_height_matrix(seq, tstruct, timepoint=5.0, use_layer_base=true)

# Create 3D mask from heights
mask_3d = create_co2_mask_3d_from_heights(height_matrix, layer, domain)

# Use mask to set CO2 saturation
co2_saturation = zeros(size(lithology))
co2_saturation[mask_3d] .= 0.6
```
"""
function create_co2_mask_3d_from_heights(
    height_matrix::Matrix{Float64},
    layer::Layer,
    domain::Domain3D
)
    nx, ny, nz = domain.nx, domain.ny, domain.nz
    dz = domain.dz
    depth_max = domain.depth_max

    # Get the layer topography (base elevation) - may be padded
    topography_2d = layer.trap_structure.topography
    pad = layer.boundary_padding

    # Validate dimensions accounting for padding
    expected_nx = nx + 2 * pad
    expected_ny = ny + 2 * pad
    @assert size(height_matrix) == (expected_nx, expected_ny) "height_matrix dimensions ($(size(height_matrix))) must match padded domain ($expected_nx, $expected_ny) for pad=$pad"

    # Initialize 3D mask (unpadded domain size)
    mask_3d = falses(nx, ny, nz)

    # Create depth array for all k indices (k=1 is deepest, k=nz is shallowest)
    depths = depth_max .- (0.5:nz) .* dz

    # For each horizontal location in the unpadded domain
    for i in 1:nx
        for j in 1:ny
            # Map to padded coordinates
            i_padded = i + pad
            j_padded = j + pad

            co2_height = height_matrix[i_padded, j_padded]

            # Skip if no CO2 at this location
            if co2_height <= 0.0
                continue
            end

            # Determine the base elevation for this location (in padded coordinates)
            base_elevation = topography_2d[i_padded, j_padded]

            # Top elevation of CO2 column at this location
            top_elevation = base_elevation + co2_height

            # Find all vertical cells within the CO2 column
            # Cell is filled if its depth is between base and top elevation
            for k in 1:nz
                cell_depth = depths[k]

                # Cell is filled if it's between base and top of CO2 column
                if base_elevation <= cell_depth <= top_elevation
                    mask_3d[i, j, k] = true
                end
            end
        end
    end

    return mask_3d
end

"""
Remove boundary padding from a 2D matrix (e.g., height matrix or topography).

This is useful for visualization and analysis after simulation with closed boundary conditions.

Parameters:
- matrix: 2D array with padding
- pad_width: Number of cells padded on each side

Returns:
- Unpadded matrix with boundary cells removed
"""
function remove_boundary_padding(matrix::Matrix{T}, pad_width::Int) where T
    if pad_width == 0
        return matrix
    end

    nx, ny = size(matrix)
    return matrix[pad_width+1:nx-pad_width, pad_width+1:ny-pad_width]
end

"""
Remove boundary padding from a vector of 2D matrices.
"""
function remove_boundary_padding(matrices::Vector{Matrix{T}}, pad_width::Int) where T
    if pad_width == 0
        return matrices
    end

    return [remove_boundary_padding(m, pad_width) for m in matrices]
end