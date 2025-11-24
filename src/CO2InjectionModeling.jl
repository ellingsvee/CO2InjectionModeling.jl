module CO2InjectionModeling

# Load subfiles
include("types.jl")
include("io.jl")
include("layer_analysis.jl")
include("volume_conversion.jl")
include("utils.jl")
include("height_tracking.jl")
include("leakage.jl")

include("visualization.jl")

include("fill_layer.jl")
include("fill_layers.jl")
include("analysis.jl")


export load_sleipner_topography, SleipnerTopography, reconstruct_3d_lithology, plot_cross_section, analyze_base_surfaces, create_trap_mask_3d, get_all_trap_masks_for_layer, Layer, Domain3D, create_domain_from_topography, footprint_to_xy, get_trap_centroid, linear_to_cartesian, animate_trap_filling, animate_trap_filling_birdseye, compute_co2_height_matrix, create_co2_mask_3d_from_heights, SimulationSummary, generate_simulation_summary

include("CO2RInterface.jl")

using .CO2RInterface
export setup_simulator, run_simulation

end # module