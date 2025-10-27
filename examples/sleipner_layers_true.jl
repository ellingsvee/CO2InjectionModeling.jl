# ============================================================================
# Sleipner CO2 Injection Simulation with Multi-Layer Leakage
# ============================================================================

using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using CairoMakie
# using GLMakie
using ColorSchemes
using NPZ
using Random
import Interpolations

# ============================================================================
# 1. Load Data
# ============================================================================

datapath = "../sleipner/sleipner_topographies.npz"
caprock_topography = NPZ.npzread(datapath)["topographies"][:, :, :]

# ============================================================================
# 2. Analyze Layer Trap Structures
# ============================================================================

tstructs = analyze_layers(caprock_topography)

# ============================================================================
# 3. Setup Injection Parameters
# ============================================================================

# Select random injection location within a trap in the top layer
top_tstruct = tstructs[1]
footprints = top_tstruct.footprints
trap_idx = rand(1:length(footprints))
footprint = footprints[trap_idx]
injection_location = rand(footprint)

# Injection rate
injection_rate = 1.0

# Compute leakage heights (half of max spillpoint height for each layer)
function compute_leakage_heights(tstructs)
    leakage_heights = Float64[]
    for i in 1:length(tstructs)-1
        tstruct = tstructs[i]
        spillpoints = tstruct.spillpoints
        max_leakage_height = -Inf

        for (j, spillpoint) in enumerate(spillpoints)
            trap_bottom = minimum(tstruct.topography[tstruct.footprints[j]])
            leakage_height = spillpoint.elevation - trap_bottom
            max_leakage_height = max(max_leakage_height, leakage_height)
        end

        push!(leakage_heights, max_leakage_height / 1.5)
    end
    return leakage_heights
end

leakage_heights = compute_leakage_heights(tstructs)


# ============================================================================
# 4. Run Simulation
# ============================================================================

seqs, times_at_leakage, leakage_locations = run_injection(
    tstructs, injection_location, injection_rate, leakage_heights
)

# Filter out nothing values from leakage times
times_at_leakage = Float64[x for x in times_at_leakage if x !== nothing]

# ============================================================================
# 5. Generate Visualization Data
# ============================================================================

# Color mapping for CO2 states
cmap = Dict(
    :blue => 2,      # Filled
    :green => 4,     # (unused)
    :red => 6,       # River
    :orange => 8,    # Trap
    :lilac => 10,    # (unused)
    :bright => 11    # (unused)
)

start_time = 0.1

# Determine how many layers actually have CO2 (filled sequences)
# CO2 reaches: layer 1 (injection) + layers where leakage occurred
n_filled_layers = length(times_at_leakage) + 1

# End time is the last timestamp of the last filled layer
end_time = seqs[n_filled_layers][end].timestamp

# Generate animation frames
times_for_animation = range(start_time, end_time, length=40)
if !isempty(times_at_leakage)
    times_for_animation = sort(collect(union(times_for_animation, times_at_leakage)))
else
    times_for_animation = collect(times_for_animation)
end

# Interpolate CO2 distribution for animation
texs_for_animation = []
for i in 1:length(tstructs)
    tstruct = tstructs[i]
    seq = seqs[i]

    # Check if this layer has any CO2 (is it within filled layers?)
    if i > n_filled_layers
        # This layer never receives CO2 - create empty texture for all times
        empty_state = ones(Int, size(tstruct.topography))  # 1 = empty state
        tex = [empty_state for _ in times_for_animation]
        push!(texs_for_animation, tex)
        continue
    end

    seq_start = seq[1].timestamp
    seq_end = seq[end].timestamp

    # Filter times based on when this layer is active
    if i <= length(times_at_leakage)
        # This layer has leakage - only animate until leakage time
        times_filtered = [t for t in times_for_animation if t <= times_at_leakage[i]]
    else
        # This is the last filled layer - animate through end
        times_filtered = times_for_animation
    end

    # Ensure times are within sequence bounds, unique, and sorted
    times_filtered = unique(sort(filter(t -> t >= seq_start && t <= seq_end, times_filtered)))
    times_filtered = length(times_filtered) >= 2 ? times_filtered : [seq_start, seq_end]

    tex, = interpolate_timeseries(tstruct, seq, times_filtered;
        filled_color=cmap[:blue],
        trap_color=cmap[:orange],
        river_color=cmap[:red])
    push!(texs_for_animation, tex)
end

# ============================================================================
# 6. Create Animation
# ============================================================================

n_layers = length(texs_for_animation)
n_cols = 3
n_rows = ceil(Int, n_layers / n_cols)

# Initialize observables for each layer
current_tex = [Observable(reverse(transpose(texs_for_animation[i][1]), dims=1)) for i in 1:n_layers]

# Create figure with better proportions
fig = Figure(size=(1600, 1300), figure_padding=20)

# Title section at the top
Label(fig[1, 1:n_cols+1], "CO2 Migration Through Sleipner Layers",
    fontsize=26, halign=:center, font=:bold, tellwidth=false)

# Time label below title
time_label = Label(fig[2, 1:n_cols+1], "Time: $(round(times_for_animation[1], digits=2)) years",
    fontsize=20, halign=:center, tellwidth=false)

# Create grid of subplots (starting at row 3)
axes = []
co2_masks = []
topo_hmaps = []

# Get elevation range across all layers for consistent colorbar
all_elevations = [caprock_topography[i, :, :] for i in 1:n_layers]
elev_min = minimum(minimum, all_elevations)
elev_max = maximum(maximum, all_elevations)

for i in 1:n_layers
    row = div(i - 1, n_cols) + 3  # Start at row 3 (after title and time)
    col = mod(i - 1, n_cols) + 1

    ax = Axis(fig[row, col];
        aspect=DataAspect(),
        title="Layer $i",
        titlesize=18,
        titlegap=8)

    # Background: topography
    topo = caprock_topography[i, :, :]
    topo_reversed = reverse(transpose(topo), dims=1)
    topo_hm = heatmap!(ax, topo_reversed;
        colormap=:terrain,
        colorrange=(elev_min, elev_max))

    # Overlay: CO2 distribution with contrasting color (use NaN for transparency)
    co2_mask = Observable(fill(NaN, size(current_tex[i][])))
    for idx in CartesianIndices(current_tex[i][])
        val = current_tex[i][][idx]
        co2_mask[][idx] = val > 1 ? val : NaN
    end

    heatmap!(ax, co2_mask;
        colormap=:hot,  # Contrasting warm colors for CO2
        colorrange=(1, 12),
        alpha=1.0,
        nan_color=:transparent)

    hidedecorations!(ax)
    push!(axes, ax)
    push!(co2_masks, co2_mask)
    push!(topo_hmaps, topo_hm)
end

# Colorbar for topography elevation
Colorbar(fig[3:n_rows+2, n_cols+1],
    topo_hmaps[1],
    label="Elevation (m)",
    labelsize=16,
    ticklabelsize=14,
    width=20)

# Animation loop
maxlength = maximum(length.(texs_for_animation))

record(fig, "../media/animation_temp.gif", 1:maxlength; framerate=5) do t
    for i in 1:n_layers
        idx = min(t, length(texs_for_animation[i]))
        new_tex = reverse(transpose(texs_for_animation[i][idx]), dims=1)
        current_tex[i][] = new_tex

        # Update CO2 mask
        for cidx in CartesianIndices(new_tex)
            val = new_tex[cidx]
            co2_masks[i][][cidx] = val > 1 ? val : NaN
        end
        notify(co2_masks[i])
    end

    # Update time label
    if t <= length(times_for_animation)
        time_label.text[] = "Time: $(round(times_for_animation[t], digits=2)) years"
    end
end

println("Animation saved to ../media/animation.gif")
