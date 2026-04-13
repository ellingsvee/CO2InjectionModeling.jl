module CO2BatchFill

using Statistics
using Printf

include("structs.jl")
include("topography.jl")
include("layer_analysis.jl")
include("unit_conversion.jl")
include("utils.jl")
include("leakage.jl")
include("fill_layer.jl")
include("fill_layers.jl")
include("analysis.jl")
include("r_interface.jl")
include("visualization.jl")
include("artifacts.jl")
include("export.jl")

end
