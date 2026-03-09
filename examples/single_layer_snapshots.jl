using CO2BatchFill
using SurfaceWaterIntegratedModeling
using SurfaceWaterIntegratedModeling: _compute_z_vol_tables
using GaussianRandomFields
using CairoMakie
using LaTeXStrings
using Random

Random.seed!(101)

# Grid / domain settings
const NX, NY = 100, 100
const LENGTH_X = 1000.0
const LENGTH_Y = 1000.0
const DX, DY = LENGTH_X / NX, LENGTH_Y / NY
const LAYER_THICK = 10.0     # m (vertical thickness of each sand layer)
const PAD_WIDTH = 2
const N_LAYERS = 3
const RESIDUAL_TRAPPING = 0.4

# Generate GRF surfaces at increasing depth
cov = CovarianceFunction(2, Matern(200, 2, σ=3.0))
pts = range(0.0, stop=(NX - 1) * DX, step=DX)

function sample_surface(base_depth)
    sample(GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)) .+ base_depth
end

surf_L1 = sample_surface(950.0)   # deepest  (injection layer)
surf_L2 = sample_surface(900.0)
surf_L3 = sample_surface(850.0)

# Build topography — sand_layers ordered shallowest-first, deepest-last
sand_layers = [
    Dict{String,Any}("name" => "Storage layer 3", "top" => surf_L3, "base" => surf_L3 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 2", "top" => surf_L2, "base" => surf_L2 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 1", "top" => surf_L1, "base" => surf_L1 .+ LAYER_THICK),
]
topo = GenericTopography(sand_layers, NX, NY, DX, DY, minimum(surf_L3), maximum(surf_L1) + LAYER_THICK)
domain = create_domain(topo, 0.1)

boundary_condition = :closed
layers = analyze_base_surfaces(topo; boundary_condition, pad_width=PAD_WIDTH)

# Reservoir properties
rp = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, 25_000.0, 5.0)
rp_caprock = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, Inf, 5.0)

println("Leakage height threshold: $(round(rp.leakage_height, digits=2)) m")

# Injection — single central well in layer 1 (deepest) only
TOTAL_RATE = 80_000.0
INJECTION_END = 10.0

rate_L1 = create_injection_rate(layers, (div(NX, 2), div(NY, 2)), TOTAL_RATE)

injection_events = [
    [InjectionEvent(0.0, rate_L1), InjectionEvent(INJECTION_END, create_injection_rate(layers, (1, 1), 0.0))],
    create_no_injection(layers),
    create_no_injection(layers),
]

# Run multi-layer simulation
println("\nRunning multi-layer fill simulation...")
seqs, leakage_states, weather_events_per_layer = fill_layers(
    layers, domain, [rp, rp, rp_caprock], injection_events; verbose=false)

# Choose which layer to plot
LAYER_IDX = 1  # Change this to 1, 2, or 3

# Find the time of first leakage from this layer
leakage_times = leakage_states[LAYER_IDX].leakage_start_time
finite_times = filter(isfinite, leakage_times)
if isempty(finite_times)
    error("No leakage occurs from layer $LAYER_IDX — try a different layer or lower the capillary pressure.")
end
t_leakage = minimum(finite_times)
println("\nFirst leakage from layer $LAYER_IDX at t = $(round(t_leakage, digits=4))")

# Three snapshot times: early plume, mid-fill, and at first leakage
t1 = t_leakage / 3
t2 = 2 * t_leakage / 3
t3 = t_leakage
timepoints = [t1, t2, t3]

println("Snapshot times: $(round.(timepoints, digits=3))")

multi_snaps = generate_multi_layer_snapshots(
    layers, seqs, leakage_states, weather_events_per_layer, timepoints)

# Build the 3-panel figure
update_theme!(
    merge(
        theme_latexfonts(),
        Theme(fontsize=25)
    )
)

max_h = ceil(rp.leakage_height)
layer = layers[LAYER_IDX]
tstruct = layer.trap_structure
pad = layer.boundary_condition == :closed ? PAD_WIDTH : 0
nx_padded, ny_padded = size(tstruct.topography)
nx = nx_padded - 2 * pad
ny = ny_padded - 2 * pad
z_vol_tables = _compute_z_vol_tables(tstruct)

x_coords = range(0.0, nx * DX, length=nx)
y_coords = range(0.0, ny * DY, length=ny)
topo_unpadded = tstruct.topography[pad+1:end-pad, pad+1:end-pad]

SHOW_EXTENTS = false  # true = binary extents, false = CO2 column heights

n_panels = 3
fig = Figure(size=(500 * n_panels, 500))
colgap!(fig.layout, 10)

for (i, snap) in enumerate(multi_snaps)
    lsnap = snap.layers[LAYER_IDX]

    h = height_map(tstruct, z_vol_tables, lsnap.trap_volumes;
        leaking=lsnap.trap_leaking, leakage_heights=lsnap.trap_leakage_height)
    h_unpadded = h[pad+1:end-pad, pad+1:end-pad]

    ax = Axis(fig[1, i];
        xlabel="x",
        ylabel=(i == 1 ? "y" : ""),
        title="t = $(round(timepoints[i], digits=2)) years",
        aspect=DataAspect(),
    )

    if SHOW_EXTENTS
        extent_array = Float64.(h_unpadded .> 0.0)
        heatmap!(ax, x_coords, y_coords, extent_array;
            colormap=[(:white, 0.0), (:dodgerblue, 0.5)],
            colorrange=(0.0, 1.0))
    else
        heatmap!(ax, x_coords, y_coords, h_unpadded;
            colormap=:Blues, colorrange=(0.0, max_h))
    end

    # Topography contours
    contour_levels = collect(range(extrema(topo_unpadded)...; length=20))
    major_levels = contour_levels[2:5:end]
    minor_levels = setdiff(contour_levels, major_levels)

    contour!(ax, x_coords, y_coords, topo_unpadded;
        levels=minor_levels, color=(:black, 1.0), linewidth=0.8, labels=false)
    contour!(ax, x_coords, y_coords, topo_unpadded;
        levels=major_levels, color=(:black, 1.0), linewidth=2.0, labels=true, labelsize=12)

    # Injection location
    scatter!(ax, [div(NX, 2) * DX], [div(NY, 2) * DY];
        color=:black, marker=:xcross, markersize=35)
end

if !SHOW_EXTENTS
    Colorbar(fig[1, 4]; colormap=:Blues, colorrange=(0.0, max_h),
        label=L"CO$_2$ column height", size=25)
end

outfile = joinpath(@__DIR__, "single_layer_snapshots_L$(LAYER_IDX).svg")
save(outfile, fig)
println("\nSaved: $outfile")
