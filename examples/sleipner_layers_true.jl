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


# Visualize all layers together
fig = Figure(size=(1200, 800))
ax = Axis3(fig[1, 1], xlabel="X", ylabel="Y", zlabel="Depth (m)", aspect=:data)

# Get a colormap for the different layers
colors = ColorSchemes.viridis[range(0, 1, length=size(caprock_topography, 1))]

for layer in 1:size(caprock_topography, 1)
    topo = caprock_topography[layer, :, :]
    # surface!(ax, topo; colormap=:viridis, alpha=0.7)
    surface!(ax, topo; color=colors[layer], alpha=0.7)
end

ax.zreversed = true

fig


# Analyze the layers to get the trap structures
tstructs = analyze_layers(caprock_topography)

# Starting the injection is the center of the top layer
# nx, ny = size(caprock_topography[1])
nx, ny = size(caprock_topography, 2), size(caprock_topography, 3)
center_ij = (div(nx, 2), div(ny, 2))
injection_location = two_dim_index_to_linear_index(caprock_topography[1, :, :], center_ij)

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
        push!(leakage_heights, max_leakage_height / 2.0) # Set leakage height to half the max leakage heighth
    end
    return leakage_heights
end

leakage_heights = compute_leakage_heights(tstructs)

# Run the injection simulation with leakage between the layers
seqs, times_at_leakage, leakage_locations = run_injection(tstructs, injection_location, injection_rate, leakage_heights)



# Remove any nothing values from tpoints
times_at_leakage = Float64[x for x in times_at_leakage if x !== nothing]

# Visualization
cmap = Dict(:blue => 2, :green => 4, :red => 6, :orange => 8,
    :lilac => 10, :bright => 11)

start_time = 0.1
end_time = seqs[end][end].timestamp

texs_at_leakage = []
for i in 1:length(tstructs)
    println("Processing layer $i for visualization...")
    tstruct = tstructs[i]
    seq = seqs[i]


    if i == length(tstructs)
        times = [start_time; times_at_leakage[1:i-1]; end_time] # Add end_time for the last layer
    else
        # Get texs at the times up to the leakage time of this specific layer
        times = [start_time; times_at_leakage[1:i]...] # Avoid exact start_time to prevent interpolation issues
    end

    tex, = interpolate_timeseries(tstruct, seq, times;
        filled_color=cmap[:blue],
        trap_color=cmap[:orange],
        river_color=cmap[:red])
    push!(texs_at_leakage, tex)
end


# For the animation we use 100 frames between the start and end time
# Also add the times at leakage to ensure they are included
times_for_animation = range(start_time, end_time, length=40)
times_for_animation = sort(collect(union(times_for_animation, times_at_leakage)))

texs_for_animation = []
for i in 1:length(tstructs)
    println("Processing layer $i for animation...")
    tstruct = tstructs[i]
    seq = seqs[i]

    times_for_animation_filtered = [t for t in times_for_animation if t <= (i < length(tstructs) ? times_at_leakage[i] : end_time)]
    tex, = interpolate_timeseries(tstruct, seq, times_for_animation_filtered;
        filled_color=cmap[:blue],
        trap_color=cmap[:orange],
        river_color=cmap[:red])
    push!(texs_for_animation, tex)
end


n_layers = length(texs_for_animation)

# Create observables for each layer
current_tex = [Observable(reverse(transpose(texs_for_animation[i][1]), dims=1)) for i in 1:n_layers]

# Build figure and axes
fig = Figure(resolution=(1000, 800))
axes = [Axis(fig[1, i]; aspect=DataAspect()) for i in 1:n_layers]

# Link images to observables
for i in 1:n_layers
    image!(axes[i], current_tex[i]; colormap=ColorSchemes.viridis, colorrange=(1, 12))
end

# Add a label for the time
time_label = Label(fig[0, :], "Time: $(round(times_for_animation[1], digits=2))", fontsize=30, halign=:center)

# Animate
maxlength = maximum(length.(texs_for_animation))
record(fig, "media/animation.gif", 1:maxlength; framerate=5) do t
    for i in 1:n_layers
        idx = min(t, length(texs_for_animation[i]))
        current_tex[i][] = reverse(transpose(texs_for_animation[i][idx]), dims=1)
    end
    time_label.text[] = "Time: $(round(times_for_animation[t], digits=2))"
end