using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using Random
using Test

Random.seed!(1)

# ---------------------------------------------------------------------------
# Helpers: create interpretable depth surfaces for testing
# ---------------------------------------------------------------------------

"""
Single Gaussian dome on a flat background.
Creates one trap at the dome center (structural high = depth minimum).
"""
function make_test_surface(nx, ny, dx, dy, base_depth;
    dome_center=(0.5, 0.5),
    dome_sigma=150.0, dome_amp=4.0)
    xs = [(i - 0.5) * dx for i in 1:nx]
    ys = [(j - 0.5) * dy for j in 1:ny]
    cx = dome_center[1] * nx * dx
    cy = dome_center[2] * ny * dy
    surface = zeros(nx, ny)
    for i in 1:nx, j in 1:ny
        surface[i, j] = base_depth
        r2 = (xs[i] - cx)^2 + (ys[j] - cy)^2
        surface[i, j] -= dome_amp * exp(-r2 / (2 * dome_sigma^2))
    end
    return surface
end

"""
Two Gaussian domes at different positions and amplitudes on a flat background.
Creates two distinct traps separated by a saddle — interpretable as two anticlines
connected by a structural low. With closed BC, produces a multi-trap hierarchy
where the shallower dome (amp2) may spill into the deeper one (amp1).
"""
function make_two_dome_surface(nx, ny, dx, dy, base_depth;
    center1=(0.3, 0.5), amp1=10.0, sigma1=100.0,
    center2=(0.7, 0.5), amp2=6.0, sigma2=100.0)
    xs = [(i - 0.5) * dx for i in 1:nx]
    ys = [(j - 0.5) * dy for j in 1:ny]
    cx1, cy1 = center1[1] * nx * dx, center1[2] * ny * dy
    cx2, cy2 = center2[1] * nx * dx, center2[2] * ny * dy
    surface = zeros(nx, ny)
    for i in 1:nx, j in 1:ny
        surface[i, j] = base_depth
        r2_1 = (xs[i] - cx1)^2 + (ys[j] - cy1)^2
        r2_2 = (xs[i] - cx2)^2 + (ys[j] - cy2)^2
        surface[i, j] -= amp1 * exp(-r2_1 / (2 * sigma1^2))
        surface[i, j] -= amp2 * exp(-r2_2 / (2 * sigma2^2))
    end
    return surface
end

# # ---------------------------------------------------------------------------
# # Test grid parameters 
# # ---------------------------------------------------------------------------
const TEST_NX, TEST_NY = 20, 20
const TEST_LENGTH_X, TEST_LENGTH_Y = 1000.0, 1000.0
const TEST_DX, TEST_DY = TEST_LENGTH_X / TEST_NX, TEST_LENGTH_Y / TEST_NY

# # Setting up the test domain
# cov = CovarianceFunction(2, Matern(100, 2, σ=1.0))
# pts_x = range(0, stop=TEST_NX * TEST_DX, step=TEST_DX)
# pts_y = range(0, stop=TEST_NX * TEST_DX, step=TEST_DY)


# # ---------------------------------------------------------------------------
# # Shared setup at module level so test files can access it.
# # ---------------------------------------------------------------------------

# const _s1 = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y)) .+ 900.0
# const _s2 = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y)) .+ 850.0


const _s1 = make_test_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 900.0)
const _s2 = make_test_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 850.0)
const _thick = 10.0


# Shallowest first, deepest last (analyze_base_surfaces reverses)
const _sand_layers = [
    Dict{String,Any}("name" => "L2", "top" => _s2, "base" => _s2 .+ _thick),
    Dict{String,Any}("name" => "L1", "top" => _s1, "base" => _s1 .+ _thick),
]

const _topo = GenericTopography(_sand_layers, TEST_NX, TEST_NY, TEST_DX, TEST_DY,
    minimum(_s2), maximum(_s1 .+ _thick))
const domain = create_domain(_topo, 1.0)
const layers = analyze_base_surfaces(_topo; boundary_condition=:closed)
const rp_leaky = ReservoirProperties(0.3, 0.2, 0.1, 12_000.0, 5.0)
const rp_sealed = ReservoirProperties(0.3, 0.2, 0.1, Inf, 5.0)
# rp_quick: very small threshold (~0.18m) to guarantee leakage in tests without large injection
const rp_quick = ReservoirProperties(0.3, 0.2, 0.1, 1000.0, 5.0)

# Set up injection events
injection_events = Vector{Vector{InjectionEvent}}()
for (i, layer) in enumerate(layers)
    topo_size = size(layer.trap_structure.topography)
    pad = layer.boundary_condition == :closed ? 1 : 0
    if i == 1
        rate = zeros(topo_size)
        rate[div(TEST_NX, 2)+pad, div(TEST_NY, 2)+pad] = 2_000.0
        push!(injection_events, [
            InjectionEvent(0.0, rate),
            InjectionEvent(10.0, zeros(topo_size)),
        ])
    else
        push!(injection_events, [InjectionEvent(0.0, zeros(topo_size))])
    end
end


# ---------------------------------------------------------------------------
# Two-dome scenario: two separate structural highs on each sand layer.
# Interpretable as two anticlines separated by a saddle.
# With closed BC, this produces a multi-trap hierarchy (2+ distinct traps).
# ---------------------------------------------------------------------------

const _td1 = make_two_dome_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 900.0)
const _td2 = make_two_dome_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 850.0)

const _td_layers = [
    Dict{String,Any}("name" => "L2", "top" => _td2, "base" => _td2 .+ _thick),
    Dict{String,Any}("name" => "L1", "top" => _td1, "base" => _td1 .+ _thick),
]
const td_topo = GenericTopography(_td_layers, TEST_NX, TEST_NY, TEST_DX, TEST_DY,
    minimum(_td2), maximum(_td1 .+ _thick))
const td_layers = analyze_base_surfaces(td_topo; boundary_condition=:closed)

# Inject at the left dome center of layer 1 (deeper trap, 30% along x)
const td_injection_events = begin
    evs = Vector{Vector{InjectionEvent}}()
    for (i, layer) in enumerate(td_layers)
        topo_size = size(layer.trap_structure.topography)
        pad = layer.boundary_condition == :closed ? 1 : 0
        if i == 1
            rate = zeros(topo_size)
            # Left dome center is at ~30% of the grid in x, 50% in y
            ix = round(Int, 0.3 * TEST_NX) + pad
            iy = div(TEST_NY, 2) + pad
            rate[ix, iy] = 2_000.0
            push!(evs, [InjectionEvent(0.0, rate), InjectionEvent(10.0, zeros(topo_size))])
        else
            push!(evs, [InjectionEvent(0.0, zeros(topo_size))])
        end
    end
    evs
end

# ---------------------------------------------------------------------------

@testset "CO2BatchFill.jl" begin
    include("test_structs.jl")
    include("test_unit_conversion.jl")
    include("test_layer_analysis.jl")
    include("test_utils.jl")
    # include("test_fill_layers.jl")
    include("test_leakage.jl")
    include("test_fill_layer.jl")
    # include("test_analysis.jl")
end
