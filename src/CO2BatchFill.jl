module CO2BatchFill

# Core types (must come first as other modules depend on them)
include("structs.jl")

# Abstract topography interface
include("topography_interface.jl")

# Core functionality
include("layer_analysis.jl")
include("volume_conversion.jl")
include("utils.jl")  # Must come before leakage.jl (provides get_all_parents)
include("leakage.jl")

include("fill_layer.jl")
include("fill_layers.jl")
include("analysis.jl")
include("visualization.jl")
include("injection_optimization.jl")

end
