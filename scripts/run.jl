using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
import Interpolations

# Setup
# Choose boundary condition: :open (CO2 can flow out) or :closed (CO2 builds up at edges)
boundary_condition = :closed  # Change to :open for open boundaries

topography = load_sleipner_topography();
domain = create_domain_from_topography(topography, 1.0);
lithology = reconstruct_3d_lithology(topography, domain);
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition);

# Create CO2 saturation field and fill a trap
layer = 1;

# Find the trap that has the largest footprint
max_size = 0
max_trap = 0
for (trap_id, footprint) in enumerate(layers[layer].trap_structure.footprints)
    global max_size, max_trap
    if length(footprint) > max_size
        max_size = length(footprint)
        max_trap = trap_id
    end
end
trap = max_trap


# co2_saturation = zeros((domain.nx, domain.ny, domain.nz));
# mask = create_trap_mask_3d(layers[layer], trap, domain);  # Layer 1, Trap 1

# # Count the number of true values in the mask
# num_cells_in_trap = count(mask)
# size(mask)

# co2_saturation[mask] .= 0.6;
# maximum(co2_saturation)

# Get the topography dimensions from the trap structure (not domain!)
trap_topo = layers[layer].trap_structure.topography
xy = CartesianIndices(size(trap_topo))[layers[layer].trap_structure.footprints[trap][40]]


reservoir_props = ReservoirProperties()

injection_rate = fill(0.0, size(layers[layer].trap_structure.topography))
# injection_rate = fill(1e6, size(layers[layer].trap_structure.topography))
injection_rate[xy] = 1e6
# injection_rate[xy[1], xy[2] + 50] = 0.5e6
injection_event_bottom_layer = [InjectionEvent(0.0, injection_rate)] # Rain with intensity 1.0 starting at time 0.0
injection_event_remaining_layers = [InjectionEvent(0.0, zeros(size(layers[layer].trap_structure.topography)))]

injection_events = Vector{Vector{InjectionEvent}}(undef, length(layers))
injection_events[1] = injection_event_bottom_layer
for i in 2:length(layers)
    injection_events[i] = injection_event_remaining_layers
end



tstruct = layers[layer].trap_structure
seqs, leakage_states, snapshots = fill_layers(
        layers,
        domain,
        reservoir_props,
        injection_events;
        num_snapshots=14,
        start_time=0.0,
        end_time=15.0,
        verbose = false
)



idx = 3
println("Snapshot at time: ", snapshots[idx].timestamp)
println("Total CO2 volume: ", snapshots[idx].total_co2_volume)


# # Create animation of trap filling (cross-section)
# animate_trap_filling(snapshots, layers[layer], lithology, domain,
#                      output_file="trap_filling_animation.gif",
#                      direction=:y,
#                      slice_index=xy[2],
#                      fps=5)

# # Create bird's eye view animation
# animate_trap_filling_birdseye(snapshots, layers[layer], domain,
#                                output_file="trap_filling_birdseye.gif",
#                                fps=5)

# Cross-section animation (all layers combined)
animate_trap_filling_multilayer(snapshots, layers, lithology, domain,
    output_file="multilayer_xsection.gif", direction=:x)

# Bird's eye grid animation (3x3 subplots)
animate_trap_filling_birdseye_multilayer(snapshots, layers, domain,
    output_file="multilayer_birdseye.gif")


timestamps = [0.01, 1.2, 1.5, 2.3]

seq_timestamps = [0.0, 1.0, 2.0, 3.0, 4.0]
Newtimepoints = [seq_timestamps[argmin(abs.(seq_timestamps .- t))] for t in timestamps]