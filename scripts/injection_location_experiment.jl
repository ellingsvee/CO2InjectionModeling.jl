using CO2InjectionModeling
using CairoMakie
using Colors
using Printf

# =============================================================================
# SETUP
# =============================================================================

println("=" ^ 70)
println("INJECTION LOCATION OPTIMIZATION EXPERIMENT")
println("=" ^ 70)

# Create output directory
mkpath("experiment_results")

# Load Sleipner reservoir
println("\nLoading Sleipner topography...")
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=:closed)
reservoir_properties = generate_reservoir_properties_for_sleipner_layers()

xy, (utm_x, utm_y, depth) = load_feeder_location(topography)
injection_events = generate_sleipner_injection_events(layers, xy)

# Get true injection location
true_location, (utm_x, utm_y, _) = load_feeder_location(topography)
println("True Sleipner injection location: $true_location")

# Sleipner historical injection parameters
total_mass_mt = 12.2

simulation_time = 15.0  # years
co2_density = reservoir_properties[1].co2_density    # kg/m³

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

"""
Run simulation with given well locations, layers, and injection rate fractions.
Returns (storage_mt, leakage_mt, seqs, leakage_states)
"""
function run_injection_simulation(
    well_locations::Vector{CartesianIndex{2}},
    well_layers::Vector{Int},
    rate_fractions::Vector{Float64},  # Fraction of total rate for each well (sum to 1)
    layers, domain, reservoir_properties;
    verbose::Bool = false
)
    n_wells = length(well_locations)
    n_layers = length(layers)
    grid_size = size(layers[1].trap_structure.topography)

    # Total injection rate (uniform over time for simplicity)
    total_rate_mt_per_year = total_mass_mt / simulation_time
    total_rate_m3_per_year = total_rate_mt_per_year * 1e9 / co2_density

    # Group wells by layer
    wells_by_layer = Dict{Int, Vector{Int}}()
    for (well_idx, layer_idx) in enumerate(well_layers)
        if !haskey(wells_by_layer, layer_idx)
            wells_by_layer[layer_idx] = Int[]
        end
        push!(wells_by_layer[layer_idx], well_idx)
    end

    # Create injection events for each layer
    zero_injection = zeros(grid_size)
    injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)
    for layer_idx in 1:n_layers
        injection_events[layer_idx] = [InjectionEvent(0.0, zero_injection)]
    end

    # Add injection for layers with wells
    for (layer_idx, well_indices) in wells_by_layer
        injection_rate_matrix = zeros(grid_size)
        for well_idx in well_indices
            loc = well_locations[well_idx]
            rate = total_rate_m3_per_year * rate_fractions[well_idx]
            injection_rate_matrix[loc] = rate
        end
        injection_events[layer_idx] = [InjectionEvent(0.0, injection_rate_matrix)]
    end

    # Run simulation
    seqs, leakage_states = fill_layers(
        layers, domain, reservoir_properties, injection_events;
        verbose=verbose
    )

    # Generate snapshots
    snapshots = generate_reservoir_snapshots(
        layers, seqs, leakage_states, domain, reservoir_properties, injection_events;
        num_snapshots=2,
        start_time=0.0,
        end_time=simulation_time,
        verbose=false
    )

    final = snapshots[end]
    storage_mt = final.total_stored_m3 * co2_density / 1e9
    leakage_mt = final.total_leaked_m3 * co2_density / 1e9

    return (storage_mt, leakage_mt, seqs, leakage_states, injection_events)
end


"""
Optimize injection rate fractions for given well locations.
"""
function optimize_rate_fractions(
    well_locations::Vector{CartesianIndex{2}},
    well_layers::Vector{Int},
    layers, domain, reservoir_properties;
    n_iterations::Int = 30,
    verbose::Bool = false
)
    n_wells = length(well_locations)

    if n_wells == 1
        return [1.0], 0.0  # Only one well, gets all injection
    end

    best_fractions = fill(1.0 / n_wells, n_wells)  # Start with equal distribution
    best_storage = -Inf

    # Simple random search for rate fractions
    for iter in 1:n_iterations
        # Generate random fractions that sum to 1
        raw = rand(n_wells)
        fractions = raw / sum(raw)

        storage, _, _, _, _ = run_injection_simulation(
            well_locations, well_layers, fractions,
            layers, domain, reservoir_properties
        )

        if storage > best_storage
            best_storage = storage
            best_fractions = fractions
            verbose && println("  Iter $iter: Storage = $(round(storage, digits=3)) Mt")
        end
    end

    return best_fractions, best_storage
end


# =============================================================================
# SCENARIO 1: BASELINE - TRUE SLEIPNER LOCATION
# =============================================================================

println("\n" * "-" ^ 70)
println("SCENARIO 1: Baseline (True Sleipner Location)")
println("-" ^ 70)

baseline_storage, baseline_leakage, baseline_seqs, baseline_leakage_states, baseline_events = run_injection_simulation(
    [true_location],
    [1],  # Layer 1
    [1.0],  # 100% to single well
    layers, domain, reservoir_properties;
    verbose=true
);

println("  Location: $true_location in Layer 1")
println("  Storage: $(round(baseline_storage, digits=3)) Mt")
println("  Leakage: $(round(baseline_leakage, digits=3)) Mt")
println("  Efficiency: $(round(100 * baseline_storage / total_mass_mt, digits=1))%")


# =============================================================================
# SCENARIO 2: OPTIMIZED SINGLE WELL
# =============================================================================

println("\n" * "-" ^ 70)
println("SCENARIO 2: Optimized Single Well Location")
println("-" ^ 70)

println("Optimizing single well location...")
single_well_result = optimize_sleipner_well_locations(
    n_wells = 1,
    allowed_layers = [1, 2, 3, 4, 5, 6, 7, 8],  # All layers except caprock
    total_mass_mt = total_mass_mt,
    algorithm = :differential_evolution,
    max_evaluations = 100,
    verbose = false,
    seed = 42
);

# Run full simulation with optimized location
opt1_storage, opt1_leakage, opt1_seqs, opt1_leakage_states, opt1_events = run_injection_simulation(
    single_well_result.optimal_locations,
    single_well_result.optimal_layers,
    [1.0],
    layers, domain, reservoir_properties
);

println("  Optimal location: $(single_well_result.optimal_locations[1]) in Layer $(single_well_result.optimal_layers[1])")
println("  Storage: $(round(opt1_storage, digits=3)) Mt")
println("  Leakage: $(round(opt1_leakage, digits=3)) Mt")
println("  Efficiency: $(round(100 * opt1_storage / total_mass_mt, digits=1))%")
println("  Improvement vs baseline: $(round(opt1_storage - baseline_storage, digits=3)) Mt ($(round(100*(opt1_storage - baseline_storage)/baseline_storage, digits=1))%)")


# =============================================================================
# SCENARIO 3: OPTIMIZED 2 WELLS
# =============================================================================

println("\n" * "-" ^ 70)
println("SCENARIO 3: Optimized 2 Wells with Rate Division")
println("-" ^ 70)

println("Step 1: Optimizing 2 well locations...")
two_well_result = optimize_sleipner_well_locations(
    n_wells = 2,
    allowed_layers = [1, 2, 3, 4, 5, 6, 7, 8],
    total_mass_mt = total_mass_mt,
    algorithm = :differential_evolution,
    max_evaluations = 100,
    verbose = false,
    seed = 123
);

println("Step 2: Optimizing injection rate division...")
opt2_fractions, _ = optimize_rate_fractions(
    two_well_result.optimal_locations,
    two_well_result.optimal_layers,
    layers, domain, reservoir_properties;
    n_iterations = 40,
    verbose = false
)

# Run full simulation with optimized locations and rates
opt2_storage, opt2_leakage, opt2_seqs, opt2_leakage_states, opt2_events = run_injection_simulation(
    two_well_result.optimal_locations,
    two_well_result.optimal_layers,
    opt2_fractions,
    layers, domain, reservoir_properties
)

println("  Well 1: $(two_well_result.optimal_locations[1]) in Layer $(two_well_result.optimal_layers[1]), rate fraction: $(round(opt2_fractions[1]*100, digits=1))%")
println("  Well 2: $(two_well_result.optimal_locations[2]) in Layer $(two_well_result.optimal_layers[2]), rate fraction: $(round(opt2_fractions[2]*100, digits=1))%")
println("  Storage: $(round(opt2_storage, digits=3)) Mt")
println("  Leakage: $(round(opt2_leakage, digits=3)) Mt")
println("  Efficiency: $(round(100 * opt2_storage / total_mass_mt, digits=1))%")
println("  Improvement vs baseline: $(round(opt2_storage - baseline_storage, digits=3)) Mt ($(round(100*(opt2_storage - baseline_storage)/baseline_storage, digits=1))%)")


# =============================================================================
# SCENARIO 4: OPTIMIZED 3 WELLS
# =============================================================================

println("\n" * "-" ^ 70)
println("SCENARIO 4: Optimized 3 Wells with Rate Division")
println("-" ^ 70)

println("Step 1: Optimizing 3 well locations...")
three_well_result = optimize_sleipner_well_locations(
    n_wells = 3,
    allowed_layers = [1, 2, 3, 4, 5, 6, 7, 8],
    total_mass_mt = total_mass_mt,
    algorithm = :differential_evolution,
    max_evaluations = 80,
    verbose = false,
    seed = 456
);

println("Step 2: Optimizing injection rate division...")
opt3_fractions, _ = optimize_rate_fractions(
    three_well_result.optimal_locations,
    three_well_result.optimal_layers,
    layers, domain, reservoir_properties;
    n_iterations = 50,
    verbose = false
)

# Run full simulation with optimized locations and rates
opt3_storage, opt3_leakage, opt3_seqs, opt3_leakage_states, opt3_events = run_injection_simulation(
    three_well_result.optimal_locations,
    three_well_result.optimal_layers,
    opt3_fractions,
    layers, domain, reservoir_properties
)

println("  Well 1: $(three_well_result.optimal_locations[1]) in Layer $(three_well_result.optimal_layers[1]), rate fraction: $(round(opt3_fractions[1]*100, digits=1))%")
println("  Well 2: $(three_well_result.optimal_locations[2]) in Layer $(three_well_result.optimal_layers[2]), rate fraction: $(round(opt3_fractions[2]*100, digits=1))%")
println("  Well 3: $(three_well_result.optimal_locations[3]) in Layer $(three_well_result.optimal_layers[3]), rate fraction: $(round(opt3_fractions[3]*100, digits=1))%")
println("  Storage: $(round(opt3_storage, digits=3)) Mt")
println("  Leakage: $(round(opt3_leakage, digits=3)) Mt")
println("  Efficiency: $(round(100 * opt3_storage / total_mass_mt, digits=1))%")
println("  Improvement vs baseline: $(round(opt3_storage - baseline_storage, digits=3)) Mt ($(round(100*(opt3_storage - baseline_storage)/baseline_storage, digits=1))%)")


# =============================================================================
# SUMMARY COMPARISON
# =============================================================================

println("\n" * "=" ^ 70)
println("SUMMARY COMPARISON")
println("=" ^ 70)

scenarios = ["Baseline (1 well)", "Optimized 1 well", "Optimized 2 wells"]
storages = [baseline_storage, opt1_storage, opt2_storage]
leakages = [baseline_leakage, opt1_leakage, opt2_leakage]
efficiencies = storages ./ total_mass_mt .* 100

println("\n" * @sprintf("%-20s %12s %12s %12s", "Scenario", "Storage(Mt)", "Leakage(Mt)", "Efficiency"))
println("-" ^ 60)
for (name, stor, leak, eff) in zip(scenarios, storages, leakages, efficiencies)
    println(@sprintf("%-20s %12.3f %12.3f %11.1f%%", name, stor, leak, eff))
end
println("-" ^ 60)


# =============================================================================
# VISUALIZATION 1: Comparison Bar Chart
# =============================================================================

println("\n" * "-" ^ 70)
println("Generating comparison visualizations...")
println("-" ^ 70)

fig_comparison = Figure(size=(1000, 600))

# Storage comparison
ax1 = Axis(fig_comparison[1, 1],
    xlabel = "Scenario",
    ylabel = "CO₂ Storage (Mt)",
    title = "Storage Comparison",
    xticks = (1:3, scenarios),
    xticklabelrotation = 0.3
)

colors_storage = [:gray, :blue, :green]
barplot!(ax1, 1:3, storages, color=colors_storage, strokewidth=1, strokecolor=:black)

# Add value labels on bars
for (i, s) in enumerate(storages)
    text!(ax1, i, s + 0.1, text="$(round(s, digits=2))", align=(:center, :bottom), fontsize=12)
end

# Efficiency comparison
ax2 = Axis(fig_comparison[1, 2],
    xlabel = "Scenario",
    ylabel = "Storage Efficiency (%)",
    title = "Efficiency Comparison",
    xticks = (1:3, scenarios),
    xticklabelrotation = 0.3
)

barplot!(ax2, 1:3, efficiencies, color=colors_storage, strokewidth=1, strokecolor=:black)

for (i, e) in enumerate(efficiencies)
    text!(ax2, i, e + 1, text="$(round(e, digits=1))%", align=(:center, :bottom), fontsize=12)
end

# Leakage comparison
ax3 = Axis(fig_comparison[2, 1:2],
    xlabel = "Scenario",
    ylabel = "CO₂ Leakage (Mt)",
    title = "Leakage Comparison",
    xticks = (1:3, scenarios)
)

barplot!(ax3, 1:3, leakages, color=[:red for _ in 1:3], strokewidth=1, strokecolor=:black, alpha=0.7)

for (i, l) in enumerate(leakages)
    text!(ax3, i, l + 0.1, text="$(round(l, digits=2))", align=(:center, :bottom), fontsize=12)
end

save("experiment_results/scenario_comparison.png", fig_comparison)
println("Saved: experiment_results/scenario_comparison.png")


# =============================================================================
# VISUALIZATION 2: Well Locations Map
# =============================================================================

fig_locations = Figure(size=(1200, 500))

# Get topography for background
topo = layers[1].trap_structure.topography
nx, ny = size(topo)

# Plot for each optimized scenario
all_locations = [
    ([true_location], [1], "Baseline"),
    (single_well_result.optimal_locations, single_well_result.optimal_layers, "Opt. 1 well"),
    (two_well_result.optimal_locations, two_well_result.optimal_layers, "Opt. 2 wells"),
    # (three_well_result.optimal_locations, three_well_result.optimal_layers, "Opt. 3 wells")
]

for (idx, (locs, lyrs, title)) in enumerate(all_locations)
    ax = Axis(fig_locations[1, idx],
        title = title,
        xlabel = idx == 1 ? "X (grid)" : "",
        ylabel = idx == 1 ? "Y (grid)" : "",
        aspect = DataAspect()
    )

    heatmap!(ax, topo, colormap=:viridis, colorrange=(minimum(topo), maximum(topo)))

    # Plot wells
    well_colors = [:red, :orange, :yellow, :cyan, :magenta]
    for (i, (loc, layer)) in enumerate(zip(locs, lyrs))
        scatter!(ax, [loc[1]], [loc[2]],
                color=well_colors[mod1(i, length(well_colors))],
                markersize=20,
                marker=:star5,
                strokewidth=2,
                strokecolor=:white)
        text!(ax, loc[1] + 2, loc[2] + 2, text="L$layer", fontsize=10, color=:white)
    end
end

save("experiment_results/well_locations_comparison.png", fig_locations)
println("Saved: experiment_results/well_locations_comparison.png")


# =============================================================================
# VISUALIZATION 3: CO2 Distribution Bird's Eye Views
# =============================================================================

println("\nGenerating CO₂ distribution plots...")

# Plot for baseline
println("  Baseline...")
plot_final_co2_distribution(
    layers,
    baseline_seqs,
    domain;
    output_file="experiment_results/co2_distribution_baseline.png",
    time=simulation_time,
    co2_color=colorant"#2166AC",  # Blue
    co2_alpha=1.0,
    show_contours=true,
    contour_levels=10,
    contour_color=:gray50,
    contour_linewidth=0.5,
    coords_in_km=true,
    transpose_layout=true,
    fontsize_layer_title=10,
    fontsize_labels=10,
    fontsize_ticks=9,
    figure_size=(900, 550)
)

# Plot for optimized single well
println("  Optimized 1 well...")
plot_final_co2_distribution(
    layers,
    opt1_seqs,
    domain;
    output_file="experiment_results/co2_distribution_opt1well.png",
    time=simulation_time,
    co2_color=colorant"#1A9850",  # Green
    co2_alpha=1.0,
    show_contours=true,
    contour_levels=10,
    contour_color=:gray50,
    contour_linewidth=0.5,
    coords_in_km=true,
    transpose_layout=true,
    fontsize_layer_title=10,
    fontsize_labels=10,
    fontsize_ticks=9,
    figure_size=(900, 550)
)

# Plot for optimized 2 wells
println("  Optimized 2 wells...")
plot_final_co2_distribution(
    layers,
    opt2_seqs,
    domain;
    output_file="experiment_results/co2_distribution_opt2wells.png",
    time=simulation_time,
    co2_color=colorant"#D73027",  # Red-orange
    co2_alpha=1.0,
    show_contours=true,
    contour_levels=10,
    contour_color=:gray50,
    contour_linewidth=0.5,
    coords_in_km=true,
    transpose_layout=true,
    fontsize_layer_title=10,
    fontsize_labels=10,
    fontsize_ticks=9,
    figure_size=(900, 550)
)

# Plot for optimized 3 wells
println("  Optimized 3 wells...")
plot_final_co2_distribution(
    layers,
    opt3_seqs,
    domain;
    output_file="experiment_results/co2_distribution_opt3wells.png",
    time=simulation_time,
    co2_color=colorant"#F46D43",  # Orange
    co2_alpha=1.0,
    show_contours=true,
    contour_levels=10,
    contour_color=:gray50,
    contour_linewidth=0.5,
    coords_in_km=true,
    transpose_layout=true,
    fontsize_layer_title=10,
    fontsize_labels=10,
    fontsize_ticks=9,
    figure_size=(900, 550)
)


# =============================================================================
# VISUALIZATION 4: Combined Summary Figure
# =============================================================================

println("\nGenerating combined summary figure...")

fig_summary = Figure(size=(1400, 900))

# Title
Label(fig_summary[0, 1:4], "Injection Location Optimization Experiment - Summary",
      fontsize=20, font=:bold)

# Row 1: Storage and efficiency bars
ax_stor = Axis(fig_summary[1, 1:2],
    ylabel = "CO₂ Storage (Mt)",
    title = "Storage by Scenario",
    xticks = (1:3, ["Baseline", "1 Well Opt.", "2 Wells Opt."])
)
barplot!(ax_stor, 1:3, storages, color=colors_storage, strokewidth=1, strokecolor=:black)
hlines!(ax_stor, [baseline_storage], color=:gray, linestyle=:dash, linewidth=2, label="Baseline")

ax_eff = Axis(fig_summary[1, 3:4],
    ylabel = "Efficiency (%)",
    title = "Storage Efficiency",
    xticks = (1:3, ["Baseline", "1 Well Opt.", "2 Wells Opt."])
)
barplot!(ax_eff, 1:3, efficiencies, color=colors_storage, strokewidth=1, strokecolor=:black)

# Row 2: Well locations
for (idx, (locs, lyrs, title)) in enumerate(all_locations)
    ax = Axis(fig_summary[2, idx],
        title = title,
        aspect = DataAspect()
    )

    heatmap!(ax, topo, colormap=:viridis)

    well_colors = [:red, :orange, :yellow, :cyan]
    for (i, (loc, layer)) in enumerate(zip(locs, lyrs))
        scatter!(ax, [loc[1]], [loc[2]],
                color=well_colors[mod1(i, length(well_colors))],
                markersize=15,
                marker=:star5,
                strokewidth=1.5,
                strokecolor=:white)
    end

    hidedecorations!(ax)
end

# Row 3: Details table (as text)
details_text = """
Scenario Details:
• Baseline: Single well at true Sleipner location ($(true_location)) in Layer 1
• 1 Well Opt.: $(single_well_result.optimal_locations[1]) in Layer $(single_well_result.optimal_layers[1])
• 2 Wells Opt.: Wells at $(two_well_result.optimal_locations[1]) (L$(two_well_result.optimal_layers[1]), $(round(opt2_fractions[1]*100, digits=0))%) and $(two_well_result.optimal_locations[2]) (L$(two_well_result.optimal_layers[2]), $(round(opt2_fractions[2]*100, digits=0))%)
"""

Label(fig_summary[3, 1:3], details_text, fontsize=11, halign=:left, valign=:top)

save("experiment_results/experiment_summary.png", fig_summary)
println("Saved: experiment_results/experiment_summary.png")


# =============================================================================
# FINAL OUTPUT
# =============================================================================

println("\n" * "=" ^ 70)
println("EXPERIMENT COMPLETE")
println("=" ^ 70)
println("\nAll results saved to experiment_results/:")
println("  - scenario_comparison.png       : Bar charts comparing scenarios")
println("  - well_locations_comparison.png : Map of well placements")
println("  - co2_distribution_baseline.png : CO₂ distribution for baseline")
println("  - co2_distribution_opt1well.png : CO₂ distribution for 1 well optimized")
println("  - co2_distribution_opt2wells.png: CO₂ distribution for 2 wells optimized")
println("  - experiment_summary.png        : Combined summary figure")

println("\nKey findings:")
println("  Best storage improvement: $(round(maximum(storages) - baseline_storage, digits=3)) Mt")
println("  Best efficiency: $(round(maximum(efficiencies), digits=1))%")
best_idx = argmax(storages)
println("  Best scenario: $(scenarios[best_idx])")
