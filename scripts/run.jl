using CO2InjectionSimulator
using SurfaceWaterIntegratedModeling
import Interpolations

# Setup
topography = load_sleipner_topography();
domain = create_domain_from_topography(topography, 1.0);
lithology = reconstruct_3d_lithology(topography, domain);
layers = analyze_base_surfaces(topography);

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


co2_saturation = zeros(size(lithology));
mask = create_trap_mask_3d(layers[layer], trap, domain);  # Layer 1, Trap 1

# Count the number of true values in the mask
num_cells_in_trap = count(mask)
size(mask)

co2_saturation[mask] .= 0.6;
maximum(co2_saturation)

# Get the topography dimensions from the trap structure (not domain!)
trap_topo = layers[layer].trap_structure.topography
xy = CartesianIndices(size(trap_topo))[layers[layer].trap_structure.footprints[trap][40]]

# # WORKS
# fig = plot_cross_section(lithology, co2_saturation, domain;
#                         slice_index=xy[1],
#                         mode=:overlay,
#                         value_label="CO2 Saturation",
#                         value_colormap=:thermal,
#                         direction=:x)

# fig = plot_cross_section(lithology, co2_saturation, domain;
#                         slice_index=xy[2],
#                         mode=:overlay,
#                         value_label="CO2 Saturation",
#                         value_colormap=:thermal,
#                         direction=:y)



reservoir_props = ReservoirProperties()

injection_rate = fill(0.0, size(layers[layer].trap_structure.topography))
# injection_rate = fill(1e6, size(layers[layer].trap_structure.topography))
injection_rate[xy] = 0.5e6
injection_rate[xy[1], xy[2] + 50] = 0.5e6
injection_events = [InjectionEvent(0.0, injection_rate)] # Rain with intensity 1.0 starting at time 0.0


seq, snapshots, leakage_state = fill_layer_co2(layers[layer].trap_structure, domain, reservoir_props, injection_events);

# Post-process: compute leaked volumes per trap
if !isnothing(leakage_state.first_leakage_time)
    # Use a time slightly after the last leakage event to capture leaked volume
    # (The simulation may end exactly when leakage occurs, giving 0 duration)
    last_leak_time = maximum([e.timestamp for e in leakage_state.leakage_events])
    final_time = max(seq[end].timestamp, last_leak_time + 0.1)  # Add 0.1 time units after last leak

    println("\n=== Leakage Analysis ===")
    println("Snapshot end time: $(snapshots[end].timestamp)")
    println("Sequence end time: $(seq[end].timestamp)")
    println("Last leak time: $last_leak_time")
    println("Final time for leakage calculation: $final_time")

    leaked_volumes_per_trap = compute_per_trap_leaked_volumes(
        leakage_state,
        injection_events,
        reservoir_props,
        domain,
        layers[layer].trap_structure,
        final_time
    )

    # Get total or per-trap volumes
    println("\nPer-trap leaked volumes:")
    for (trap_id, vol) in sort(collect(leaked_volumes_per_trap))
        vol_physical = swim_volume_to_physical_volume(vol, reservoir_props, domain)
        println("  Trap $trap_id: $(round(vol, digits=2)) swim units = $(round(vol_physical, digits=2)) m³")
    end

    total_leaked = get_total_leaked_volume(leaked_volumes_per_trap)
    total_physical = swim_volume_to_physical_volume(total_leaked, reservoir_props, domain)
    println("\nTotal leaked: $(round(total_leaked, digits=2)) swim units = $(round(total_physical, digits=2)) m³")

    # Debug: Which region/trap contains xy?
    num_regions = SurfaceWaterIntegratedModeling.numregions(layers[layer].trap_structure)
    num_traps = SurfaceWaterIntegratedModeling.numtraps(layers[layer].trap_structure)
    linear_xy = LinearIndices(size(layers[layer].trap_structure.topography))[xy]

    println("\nDebug: num_regions=$num_regions, num_traps=$num_traps, xy=$xy, linear_xy=$linear_xy")

    # Check all traps
    found_region = nothing
    for trap_id in 1:num_traps
        if linear_xy in layers[layer].trap_structure.footprints[trap_id]
            if trap_id <= num_regions
                println("  xy is in region $trap_id")
                found_region = trap_id
            else
                println("  xy is in parent trap $trap_id")
            end
        end
    end

    if isnothing(found_region)
        println("  ERROR: xy is not in any trap footprint!")
        # Check total domain coverage
        total_covered = Set{Int}()
        for trap_id in 1:num_regions
            union!(total_covered, layers[layer].trap_structure.footprints[trap_id])
        end
        println("  Regions cover $(length(total_covered)) / $(prod(size(layers[layer].trap_structure.topography))) cells")
    end
else
    println("\n=== No leakage occurred ===")
end



# @code_warntype fill_layer_co2(layers[layer].trap_structure, domain, reservoir_props, injection);
# @time fill_layer_co2(layers[layer].trap_structure, domain, reservoir_props, injection);



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





# # For each snapshot, compute the height matrix
# i = length(snapshots)
# height_matrices = compute_co2_height_matrix(seq, layers[layer].trap_structure; timepoint = [snapshots[i].timestamp], use_layer_base = true)





    
# maximum(height_matrices[1])

# height_matrix = height_matrices[1]

# # find the index with the maximum height
# max_height = maximum(height_matrix)
# max_index = findfirst(==(max_height), height_matrix)

# height_matrix[max_index[1]-2:max_index[1]+2, max_index[2]-2:max_index[2]+2]