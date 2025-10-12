using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
import GLMakie, Images # for visualization and loading of textures
import ColorSchemes
import Graphs
import Interpolations
using Pkg.Artifacts

datapath = joinpath(datapath_testdata(), "data", "synthetic")
grid = loadgrid(joinpath(datapath, "synsurf.txt"))

cmap = Dict(:blue => 2, :green => 4, :red => 6, :orange => 8,
    :lilac => 10, :bright => 11)

## We also define some view angles that will come in handy
view1 = (GLMakie.Vec(-83, 378, 197), GLMakie.Vec(100, 114, -4.5), 0.68)
view2 = (GLMakie.Vec(90, 494, 8.6), GLMakie.Vec(100, 114, -4.5), 0.68);

## plot the grid
sf, fig, sc = plotgrid(grid, texture=fill(cmap[:bright], size(grid)),
    colormap=ColorSchemes.:Paired_12,
    colorrange=(1, 12), wireframe=true)

tstruct = spillanalysis(grid);

# weather = [WeatherEvent(0.0, 1.0)] # Rain with intensity 1.0 starting at time 0.0
injection_location = tstruct.footprints[1][1]
rain_rate_matrix = fill(0.0, size(grid))
rain_rate_matrix[injection_location] = 1.0
weather = [WeatherEvent(0.0, rain_rate_matrix)]

seq = fill_sequence(tstruct, weather);

for i = 1:length(seq)
    print("Time: ", seq[i].timestamp, ", trap fill states: ", filled_at(seq, i), '\n')
end
# amount_at(seq, 4) ## query the third event in the sequence

leakage_height = 5.0
heights_above_trap_bottom = spillpoint_heights_above_trap_bottom(tstruct)

time_at_leakage = get_time_of_leakage(seq, tstruct, injection_location, leakage_height)
if isnothing(time_at_leakage)
    println("No leakage occurred.")
end

tstates = trap_states_at_timepoints(tstruct, seq, [time_at_leakage])
water_content = [e[2] for e in tstates]
for time_ix = 1:length(tpoints)
    print("At time: ", tpoints[time_ix], ":\n")
    for trap_ix = 1:4
        content = water_content[time_ix][trap_ix]
        subtraps = Graphs.inneighbors(tstruct.agglomerations, trap_ix)
        for i in subtraps
            # add in water from subtraps
            content += water_content[time_ix][i]
        end
        print("   Trap: ", trap_ix, " contains: ", content, " units of water.\n")
    end
end