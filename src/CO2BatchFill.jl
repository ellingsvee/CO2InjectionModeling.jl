module CO2BatchFill

using Statistics
using Printf


# Core types
include("structs.jl")

# Topography
include("topography.jl")

# Core functionality
include("layer_analysis.jl")
include("unit_conversion.jl")
include("utils.jl")
include("leakage.jl")
include("fill_layer.jl")
include("fill_layers.jl")
include("analysis.jl")

# R interface
include("r_interface.jl")

# Visualization function stubs (implementations provided by extension)
include("visualization.jl")

end
