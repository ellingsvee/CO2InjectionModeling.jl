using CO2InjectionModeling
using SurfaceWaterIntegratedModeling
import GLMakie, Images # for visualization and loading of textures
import ColorSchemes
import Graphs
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

function get_dist_to_lowest_point(tstruct)
    N = numtraps(tstruct)

    dists = zeros(Float64, N)
    for i = 1:N
        spillpoint = tstruct.spillpoints[i]
        if !isnothing(spillpoint)
            lowest_point = minimum(tstruct.topography[tstruct.footprints[i]])
            dists[i] = spillpoint.elevation - lowest_point
        else
            dists[i] = 0.0
        end
    end

    return dists
end

dists = get_dist_to_lowest_point(tstruct)

# weather = [WeatherEvent(0.0, 1.0)] # Rain with intensity 1.0 starting at time 0.0
injection_location = tstruct.footprints[2][1]
rain_rate_matrix = fill(0.0, size(grid))
rain_rate_matrix[injection_location] = 1.0
weather = [WeatherEvent(0.0, rain_rate_matrix)]

seq = fill_sequence(tstruct, weather);

for i = 1:length(seq)
    print("Time: ", seq[i].timestamp, ", trap fill states: ", filled_at(seq, i), '\n')
end
# amount_at(seq, 4) ## query the third event in the sequence

leakage_height = 5.0
leakage_root_trap = get_leakage_root_trap(tstruct, injection_location, leakage_height)
children_filled_idx = get_seq_index_when_children_are_filled(tstruct, seq, 4)


out = get_time_of_leakage(tstruct, seq, leakage_height, leakage_root_trap, children_filled_idx)




tpoints = [seq[filled_idx].timestamp]
tstates = trap_states_at_timepoints(tstruct, seq, tpoints)
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


# amounts = [e.amount for e ∈ amount_at(seq, 3)] ## amounts when last computed
# dt = seq[3].timestamp .- [e.time for e ∈ amount_at(seq, 3)] ## time since last update
# amounts += dt .* inflow_at(seq, 3) # add inflow since last update

# tpoints = [0.1, 0.3, 0.42]

# # We can now compute the water content in all traps at this timepoint as
# # follows:
# tstates = trap_states_at_timepoints(tstruct, seq, tpoints)
# tex, = interpolate_timeseries(tstruct, seq, tpoints,
#     filled_color=cmap[:blue],
#     trap_color=cmap[:orange],
#     river_color=cmap[:red])

# drape_surface(sf, tex[1])

# set_camerapos(sc, view1...)
# # At time 0.1

# # drape_surface(sf, tex[2])
# # set_camerapos(sc, view1...)
# # # At time 0.3

# # drape_surface(sf, tex[3])
# # set_camerapos(sc, view1...)