using CO2BatchFill
using SurfaceWaterIntegratedModeling
using Printf

# Setup similar to debugging_single_layer.jl
boundary_condition = :closed
topography = load_sleipner_topography("examples/Sleipner/depth_surfaces")
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

# Load injection location from feeder chimney data
xy, (utm_x, utm_y, depth) = load_feeder_location(topography)
injection_events = generate_sleipner_injection_events(layers, xy)

layer_idx = 1
rprops = generate_fitted_reservoir_properties_for_sleipner_layers()

# Run fill_layer with leakage detection
seqs, leakage_states = fill_layers(
    layers,
    domain,
    rprops,
    injection_events;
    verbose=false
);

# Plotting results
using CairoMakie
using Colors
plot_final_co2_distribution(
    layers,
    seqs,
    domain;
    output_file="final_co2_distribution.svg",
    time=15.0,
    co2_color=colorant"#A49841",
    co2_alpha=1.0,
    show_contours=true,
    contour_levels=12,
    contour_color=:gray50,
    contour_linewidth=0.8,
    coords_in_km=true,  # Compact coordinate display (0-3km instead of large UTM values)
    transpose_layout=true,  # Rotate plots 90° for wider aspect
    fontsize_layer_title=12,
    fontsize_labels=12,
    fontsize_ticks=12,
    figure_size=(800, 500)
)
