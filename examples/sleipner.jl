using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
import GLMakie, Images # for visualization and loading of textures
# import CairoMakie, Images # for visualization and loading of textures
import ColorSchemes
import Graphs
import Interpolations
using Pkg.Artifacts

# datapath = joinpath(datapath_testdata(), "data", "synthetic")
# grid = loadgrid(joinpath(datapath, "synsurf.txt"))

using NPZ
datapath = "sleipner/caprock_topography.npy"
caprock_topography = NPZ.npzread(datapath)

caprock_topography = caprock_topography * 200.0
caprock_topography .-= minimum(caprock_topography) # set lowest point to zero
grid = caprock_topography


cmap = Dict(:blue => 2, :green => 4, :red => 6, :orange => 8,
    :lilac => 10, :bright => 11)

## We also define some view angles that will come in handy
view1 = (GLMakie.Vec(-83, 378, 197), GLMakie.Vec(100, 114, -4.5), 0.68)
view2 = (GLMakie.Vec(90, 494, 8.6), GLMakie.Vec(100, 114, -4.5), 0.68);
view3 = (GLMakie.Vec(-83, 378, 197), GLMakie.Vec(100, 114, -4.5), 1.2)


# view1 = (CairoMakie.Vec(-83, 378, 197), CairoMakie.Vec(100, 114, -4.5), 0.68)
# view2 = (CairoMakie.Vec(90, 494, 8.6), CairoMakie.Vec(100, 114, -4.5), 0.68);
# view3 = (CairoMakie.Vec(0, 0, 500), CairoMakie.Vec(100, 114, -4.5), 0.68);

## plot the grid
sf, fig, sc = plotgrid(grid, texture=fill(cmap[:bright], size(grid)),
    colormap=ColorSchemes.:Paired_12,
    colorrange=(1, 12), wireframe=true)


drape_surface(sf, zero(grid))

set_camerapos(sc, view3...)






tstruct = spillanalysis(grid);

# weather = [WeatherEvent(0.0, 1.0)] # Rain with intensity 1.0 starting at time 0.0
# injection_location = tstruct.footprints[1][1]
nx, ny = size(grid)
center_ij = (div(nx, 2), div(ny, 2))
# injection_location = LinearIndices(grid)[center_ij...]
injection_location = two_dim_index_to_linear_index(grid, center_ij)

# injection_location = 386888

rain_rate_matrix = fill(0.0, size(grid))
rain_rate_matrix[injection_location] = 5000.0
weather = [WeatherEvent(0.0, rain_rate_matrix)]

seq = fill_sequence(tstruct, weather);

for i = 1:length(seq)
    print("Time: ", seq[i].timestamp, ", trap fill states: ", filled_at(seq, i), '\n')
end


leakage_height = 10.0
heights_above_trap_bottom = spillpoint_heights_above_trap_bottom(tstruct)


(leakage_location, time_at_leakage) = get_loc_and_time_of_leakage(seq, tstruct, injection_location, leakage_height)
if isnothing(time_at_leakage)
    println("No leakage occurred.")
end


tpoints = [time_at_leakage * 0.25, time_at_leakage * 0.5, time_at_leakage * 0.75, time_at_leakage]
tex, = interpolate_timeseries(tstruct, seq, tpoints,
    filled_color=cmap[:blue],
    trap_color=cmap[:orange],
    river_color=cmap[:red])

drape_surface(sf, tex[5])
set_camerapos(sc, view3...)