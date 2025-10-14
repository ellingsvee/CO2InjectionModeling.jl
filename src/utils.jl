export two_dim_index_to_linear_index
export linear_index_to_two_dim_index


function two_dim_index_to_linear_index(grid, idxs)
    return LinearIndices(grid)[idxs...]
end

function linear_index_to_two_dim_index(grid, idx)
    return Tuple(CartesianIndices(grid)[idx])
end