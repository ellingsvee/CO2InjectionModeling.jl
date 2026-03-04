using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using Makie
using CairoMakie

using CO2BatchFill
using Random

Random.seed!(4)

# ---------------------------------------------------------------------------
# Grid / domain settings  (mirrors the test setup)
# ---------------------------------------------------------------------------
const NX, NY = 200, 200
const LENGTH_X = 1000.0   # m
const LENGTH_Y = 1000.0   # m
const DX, DY = LENGTH_X / NX, LENGTH_Y / NY
const LAYER_THICKNESS = 10.0     # m
const pad_width = 2

# ---------------------------------------------------------------------------
# GRF topography — same covariance model as in runtests.jl
# ---------------------------------------------------------------------------
cov = CovarianceFunction(2, Matern(100, 2, σ=1.0))
pts_x = range(0.0, stop=(NX - 1) * DX, step=DX)
pts_y = range(0.0, stop=(NY - 1) * DY, step=DY)

grf_surface = sample(
    GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y, minpadding=113)
) .+ 900.0  # mean depth ≈ 900 m

# function _make_dome(nx, ny, dx, dy, base_depth;
#     dome_center=(0.5, 0.5),
#     dome_sigma=150.0, dome_amp=4.0)
#     xs = [(i - 0.5) * dx for i in 1:nx]
#     ys = [(j - 0.5) * dy for j in 1:ny]
#     cx = dome_center[1] * nx * dx
#     cy = dome_center[2] * ny * dy
#     surface = zeros(nx, ny)
#     for i in 1:nx, j in 1:ny
#         surface[i, j] = base_depth
#         r2 = (xs[i] - cx)^2 + (ys[j] - cy)^2
#         surface[i, j] -= dome_amp * exp(-r2 / (2 * dome_sigma^2))
#     end
#     return surface
# end

# grf_surface = _make_dome(NX, NY, DX, DY, 900.0)




# ---------------------------------------------------------------------------
# Build topography and analyze trap structure
# ---------------------------------------------------------------------------
layers_raw = [
    Dict{String,Any}("name" => "Sand", "top" => grf_surface,
        "base" => grf_surface .+ LAYER_THICKNESS),
]
topo = GenericTopography(layers_raw, NX, NY, DX, DY,
    minimum(grf_surface),
    maximum(grf_surface) + LAYER_THICKNESS)
domain = create_domain(topo, 1.0)

boundary_condition = :closed
layers = analyze_base_surfaces(topo; boundary_condition, pad_width=pad_width)
layer = layers[1]   # only one layer

# ---------------------------------------------------------------------------
# Reservoir properties  (sealed caprock so CO2 accumulates in traps)
# ---------------------------------------------------------------------------
rp = ReservoirProperties(0.3, 0.2, 0.1, Inf, 5.0)
# rp = ReservoirProperties(0.3, 0.2, 0.1, 10_000.0, 5.0)
println(rp.leakage_height)

# ---------------------------------------------------------------------------
# Injection — single central well
# ---------------------------------------------------------------------------
TOTAL_RATE = 20_000.0   # m³/day
INJECTION_END_T = 15.0

pad = boundary_condition == :closed ? pad_width : 0
topo_size = (NX + 2 * pad, NY + 2 * pad)

rate = zeros(topo_size)
# rate[div(NX, 2)+pad, div(NY, 2)+pad] = TOTAL_RATE
rate[1+pad, 1+pad] = TOTAL_RATE
rate[NX-pad, NY-pad] = TOTAL_RATE

injection_events = [
    InjectionEvent(0.0, rate),
    InjectionEvent(INJECTION_END_T, zeros(topo_size)),
]

weather_events = convert_injection_event_to_weather_event(injection_events, rp, domain)

# ---------------------------------------------------------------------------
# Run simulation
# ---------------------------------------------------------------------------
println("Running fill simulation...")
seq, leakage_state = fill_sequence_with_leakage(
    layer.trap_structure, rp, weather_events; verbose=true)
println("Simulation complete — $(length(seq)) spill events.")

# ---------------------------------------------------------------------------
# Animate
# ---------------------------------------------------------------------------
t_end = maximum(se.timestamp for se in seq if isfinite(se.timestamp))

animate_layer_filling(
    layer, seq, domain;
    output_file=joinpath(@__DIR__, "single_layer_filling.gif"),
    num_frames=80,
    end_time=t_end,
    fps=15,
    max_co2_height=8.0,
    show_contours=true,
    contour_levels=12,
);
