using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using CairoMakie
using GLMakie
using ColorSchemes
import Graphs
import Interpolations
using Pkg.Artifacts

# Import and process the data
using NPZ
datapath = "../sleipner/sleipner_topographies.npz"
caprock_topography = NPZ.npzread(datapath)["topographies"][:, :, :]

#####################################################

# Visualize all layers together
fig = Figure(size=(1200, 800))
ax = Axis3(fig[1, 1], xlabel="X", ylabel="Y", zlabel="Depth (m)", aspect=:data)

for layer in 1:size(caprock_topography, 1)
    topo = caprock_topography[layer, :, :]
    surface!(ax, topo; colormap=:viridis, alpha=1)
end

ax.zreversed = true

fig

#####################################################

# 2D contour visualization
fig_contour = Figure(size=(1400, 1000))

# 2D visualization of all layers in a grid
n_layers = size(caprock_topography, 1)
n_cols = 3
n_rows = ceil(Int, n_layers / n_cols)

for layer in 1:n_layers
    topo = caprock_topography[layer, :, :]
    row = div(layer - 1, n_cols) + 1
    col = mod(layer - 1, n_cols) + 1

    ax = Axis(fig_contour[row, col],
        title="Layer $layer - Contours",
        xlabel="X",
        ylabel="Y",
        aspect=DataAspect())

    contourf!(ax, topo; colormap=:terrain, levels=15)
    contour!(ax, topo; color=:black, linewidth=0.5, levels=10)
end

fig_contour

###################################################

# Analyze the layers to get the trap structures
tstructs = analyze_layers(caprock_topography)

##################################################

layer_idx = 3
linear_idx = 2319

# Get the topography for layer 2
topo = caprock_topography[layer_idx, :, :]

# Plot footprints as contour outlines
fig_contours = Figure(size=(800, 600))
ax = Axis(fig_contours[1, 1],
    title="Layer $layer_idx - Traps as contours",
    xlabel="X",
    ylabel="Y",
    aspect=DataAspect())

# Background elevation
heatmap!(ax, topo; colormap=:terrain)

# Draw contours around each footprint
for (trap_id, footprint) in enumerate(tstruct.footprints)
    # Create binary mask for this footprint
    mask = zeros(Bool, ny, nx)
    for linear_idx in footprint
        j = div(linear_idx - 1, ny) + 1
        i = mod(linear_idx - 1, ny) + 1
        mask[i, j] = true
    end

    # Draw contour
    contour!(ax, Float64.(mask); levels=[0.5], color=:red, linewidth=2)
end

i, j = linear_index_to_two_dim_index(topo, linear_idx)

scatter!(ax, [j], [i]; color=:white, markersize=10, marker=:x, strokewidth=3)

fig_contours

###################################################

top_tstruct = tstructs[1]
footprints = top_tstruct.footprints

# Select a random loaction within one of the footprints
using Random
trap_idx = rand(1:length(footprints))
footprint = footprints[trap_idx]
injection_location = rand(footprint)



# Starting the injection is the center of the top layer
# # nx, ny = size(caprock_topography[1])
# nx, ny = size(caprock_topography, 2), size(caprock_topography, 3)
# center_ij = (div(nx, 2), div(ny, 2))
# injection_location = two_dim_index_to_linear_index(caprock_topography[1, :, :], center_ij)

# Injection rate and leakage heights
injection_rate = 1.0

function compute_leakage_heights(tstructs)
    leakage_heights = Float64[]
    for i in 1:length(tstructs)-1
        tstruct = tstructs[i]
        spillpoints = tstruct.spillpoints

        max_leakage_height = -Inf

        # Get the spillpoint with the lowest elevation
        for (i, spillpoint) in enumerate(spillpoints)

            trap_bottom = minimum(tstruct.topography[tstruct.footprints[i]])
            leakage_height = spillpoint.elevation - trap_bottom
            if leakage_height > max_leakage_height
                max_leakage_height = leakage_height
            end
        end
        push!(leakage_heights, max_leakage_height / 2) # Set leakage height to half the max leakage heighth
    end
    return leakage_heights
end

leakage_heights = compute_leakage_heights(tstructs)


##################################################

# Run the injection simulation with leakage between the layers
seqs, times_at_leakage, leakage_locations = run_injection(tstructs, injection_location, injection_rate, leakage_heights)


##################################################

# Remove any nothing values from times_at_leakage
times_at_leakage = Float64[x for x in times_at_leakage if x !== nothing]

# Visualization
cmap = Dict(:blue => 2, :green => 4, :red => 6, :orange => 8,
    :lilac => 10, :bright => 11)

start_time = 0.1

texs_at_leakage = []
for i in 1:length(tstructs)
    println("Processing layer $i for visualization...")
    tstruct = tstructs[i]
    seq = seqs[i]

    # Get the actual end time for this sequence
    seq_end_time = seq[end].timestamp

    # Build times array based on available leakage times
    if isempty(times_at_leakage)
        # No leakage occurred - just use start and end times
        times = [start_time, seq_end_time]
    elseif i == length(tstructs)
        # Last layer - include all previous leakage times plus end time
        times = [start_time; times_at_leakage; seq_end_time]
    elseif i <= length(times_at_leakage)
        # Include leakage times up to this layer
        times = [start_time; times_at_leakage[1:i]...]
    else
        # This layer is beyond where leakage occurred
        times = [start_time; times_at_leakage; seq_end_time]
    end

    # Ensure times are within sequence bounds and unique
    times = unique(sort(filter(t -> t >= seq[1].timestamp && t <= seq_end_time, times)))

    tex, = interpolate_timeseries(tstruct, seq, times;
        filled_color=cmap[:blue],
        trap_color=cmap[:orange],
        river_color=cmap[:red])
    push!(texs_at_leakage, tex)
end

# Get overall end time from the last sequence
end_time = seqs[end][end].timestamp


# For the animation we use 40 frames between the start and end time
# Also add the times at leakage to ensure they are included
times_for_animation = range(start_time, end_time, length=40)
if !isempty(times_at_leakage)
    times_for_animation = sort(collect(union(times_for_animation, times_at_leakage)))
else
    times_for_animation = collect(times_for_animation)
end

texs_for_animation = []
for i in 1:length(tstructs)
    println("Processing layer $i for animation...")
    tstruct = tstructs[i]
    seq = seqs[i]

    # Get sequence bounds
    seq_start = seq[1].timestamp
    seq_end = seq[end].timestamp

    # Filter times based on when leakage occurs (if at all)
    if isempty(times_at_leakage)
        # No leakage - use all animation times for all layers
        times_for_animation_filtered = times_for_animation
    elseif i < length(tstructs) && i <= length(times_at_leakage)
        # Filter up to this layer's leakage time
        times_for_animation_filtered = [t for t in times_for_animation if t <= times_at_leakage[i]]
    else
        # Last layer or beyond leakage - use all times
        times_for_animation_filtered = times_for_animation
    end

    # Ensure times are within sequence bounds, unique, and sorted
    times_for_animation_filtered = unique(sort(filter(t -> t >= seq_start && t <= seq_end, times_for_animation_filtered)))

    # Ensure we have at least two points
    if length(times_for_animation_filtered) < 2
        times_for_animation_filtered = [seq_start, seq_end]
    end

    tex, = interpolate_timeseries(tstruct, seq, times_for_animation_filtered;
        filled_color=cmap[:blue],
        trap_color=cmap[:orange],
        river_color=cmap[:red])
    push!(texs_for_animation, tex)
end


n_layers = length(texs_for_animation)

# Create observables for each layer
current_tex = [Observable(reverse(transpose(texs_for_animation[i][1]), dims=1)) for i in 1:n_layers]

# Grid layout parameters (3x3 for 9 layers)
n_cols = 3
n_rows = ceil(Int, n_layers / n_cols)

# Build figure with grid layout
fig = Figure(size=(1400, 1200))

# Add title
Label(fig[0, :], "CO2 Migration Through Sleipner Layers", fontsize=24,
    halign=:center, font=:bold)

# Add time label
time_label = Label(fig[1, :], "Time: $(round(times_for_animation[1], digits=2)) years",
    fontsize=20, halign=:center)

# Create axes in grid layout
axes = []
topo_plots = []

for i in 1:n_layers
    row = div(i - 1, n_cols) + 1 + 1  # +1 to account for title row
    col = mod(i - 1, n_cols) + 1

    ax = Axis(fig[row+1, col];
        aspect=DataAspect(),
        title="Layer $i",
        titlesize=16)

    # Plot topography as background with terrain colormap
    topo = caprock_topography[i, :, :]
    topo_reversed = reverse(transpose(topo), dims=1)

    # Show topography with terrain colors
    topo_hm = heatmap!(ax, topo_reversed;
        colormap=:terrain,
        alpha=1.0)

    # Create a mask for CO2 visualization
    # We'll only show CO2 where it exists (not in empty cells)
    co2_mask = Observable(similar(current_tex[i][], Float64))
    co2_mask[] .= NaN  # Start with all NaN (transparent)

    # Update mask based on CO2 state
    for idx in CartesianIndices(current_tex[i][])
        val = current_tex[i][][idx]
        if val > 1  # If not empty (1 = empty)
            co2_mask[][idx] = val
        else
            co2_mask[][idx] = NaN  # Transparent for empty cells
        end
    end

    # Overlay CO2 distribution (only where CO2 exists)
    img = heatmap!(ax, co2_mask;
        colormap=ColorSchemes.viridis,
        colorrange=(1, 12),
        alpha=0.7,
        nan_color=:transparent)

    # Hide decorations for cleaner look
    hidedecorations!(ax)

    push!(axes, ax)
    push!(topo_plots, co2_mask)
end

# Add colorbar for CO2 state
Colorbar(fig[2:n_rows+1, n_cols+1],
    limits=(1, 12),
    colormap=ColorSchemes.viridis,
    label="CO2 State",
    labelsize=14,
    ticklabelsize=12,
    ticks=(1:12, ["Empty", "Filled", "Trap", "River",
        "5", "6", "7", "8", "9", "10", "11", "12"]))

# Animate
maxlength = maximum(length.(texs_for_animation))
println("Creating animation with $maxlength frames...")

record(fig, "../media/animation_true.gif", 1:maxlength; framerate=5) do t
    for i in 1:n_layers
        idx = min(t, length(texs_for_animation[i]))
        new_tex = reverse(transpose(texs_for_animation[i][idx]), dims=1)
        current_tex[i][] = new_tex

        # Update CO2 mask to only show where CO2 exists
        for cidx in CartesianIndices(new_tex)
            val = new_tex[cidx]
            if val > 1  # If not empty
                topo_plots[i][][cidx] = val
            else
                topo_plots[i][][cidx] = NaN  # Transparent for empty
            end
        end
        notify(topo_plots[i])
    end

    # Update time label
    if t <= length(times_for_animation)
        time_label.text[] = "Time: $(round(times_for_animation[t], digits=2)) years"
    end
end

println("Animation saved to ../media/animation_true.gif")