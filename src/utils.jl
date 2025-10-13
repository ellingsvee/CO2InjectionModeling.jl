export two_dim_index_to_linear_index

function two_dim_index_to_linear_index(grid, idxs)
    return LinearIndices(grid)[idxs...]
end

