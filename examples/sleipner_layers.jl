using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
using CairoMakie, ColorSchemes
import Graphs
import Interpolations
using Pkg.Artifacts

# Import and process the data
using NPZ
datapath = "sleipner/caprock_topography.npy"
caprock_topography = NPZ.npzread(datapath)
caprock_topography = caprock_topography * 200.0
caprock_topography .-= minimum(caprock_topography) # set lowest point to zero
caprock_topography_rot_180 = reverse(caprock_topography, dims=(1, 2))

# Spesify the layers
dist_between_layers = 10.0 + maximum(caprock_topography)
layer1 = caprock_topography .+ 2 * dist_between_layers
layer2 = caprock_topography_rot_180 .+ dist_between_layers
layer3 = caprock_topography
layers = [layer1, layer2, layer3]

# Analyze the layers to get the trap structures
tstructs = analyze_layers(layers)

# Starting the injection is the center of the top layer
nx, ny = size(layers[1])
center_ij = (div(nx, 2), div(ny, 2))
injection_location = two_dim_index_to_linear_index(layers[1], center_ij)

# Injection rate and leakage heights
injection_rate = 5000.0
leakage_heights = [20.0, 15.0] # leakage heights

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