using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using Random
using Test

Random.seed!(1)

# ---------------------------------------------------------------------------
# Test grid parameters 
# ---------------------------------------------------------------------------
const TEST_NX, TEST_NY = 200, 200
const TEST_LENGTH_X, TEST_LENGTH_Y = 1000.0, 1000.0
const TEST_DX, TEST_DY = TEST_LENGTH_X / TEST_NX, TEST_LENGTH_Y / TEST_NY
const _thick = 10.0
const boundary_condition = :closed # :open or :closed

const TOTAL_INJECTION_RATE = 20_000.0
const INJECTION_START_T = 0.0
const INJECTION_END_T = 15.0

const pad_width = 2

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function create_domain_and_layers(
    top_layer::Matrix{<:Real}, bottom_layer::Matrix{<:Real}, thickness::Real,
    nx::Int, ny::Int, dx::Float64, dy::Float64,
    boundary_condition::Symbol
)

    # Shallowest first, deepest last (analyze_base_surfaces reverses)
    _layers = [
        Dict{String,Any}("name" => "L2", "top" => top_layer, "base" => top_layer .+ thickness),
        Dict{String,Any}("name" => "L1", "top" => bottom_layer, "base" => bottom_layer .+ thickness),
    ]
    topo = GenericTopography(_layers, nx, ny, dx, dy, minimum(top_layer), maximum(bottom_layer) + thickness)
    domain = create_domain(topo, 1.0)
    analyzed_layers = analyze_base_surfaces(topo; boundary_condition=boundary_condition, pad_width=pad_width)
    return topo, domain, analyzed_layers
end

# ---------------------------------------------------------------------------
# Single-dome scenario: one structural high on each sand layer.
# ----------------------------------------------------------------------------

function _make_dome(nx, ny, dx, dy, base_depth;
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

const _dome1 = _make_dome(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 850.0)
const _dome2 = _make_dome(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 900.0)

dome_topo, dome_domain, dome_layers = create_domain_and_layers(_dome1, _dome2, _thick, TEST_NX, TEST_NY, TEST_DX, TEST_DY, boundary_condition)


# ---------------------------------------------------------------------------
# Realistic scenario: Topography sampled from a GRF
# ---------------------------------------------------------------------------
cov = CovarianceFunction(2, Matern(100, 2, σ=1.0))
pts_x = range(0, stop=(TEST_NX - 1) * TEST_DX, step=TEST_DX)
pts_y = range(0, stop=(TEST_NY - 1) * TEST_DY, step=TEST_DY)

const _grf1 = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y, minpadding=113)) .+ 900.0
const _grf2 = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y, minpadding=113)) .+ 850.0

grf_topo, grf_domain, grf_layers = create_domain_and_layers(_grf1, _grf2, _thick, TEST_NX, TEST_NY, TEST_DX, TEST_DY, boundary_condition)

# ---------------------------------------------------------------------------
# Reservoir properties
# ---------------------------------------------------------------------------

const rp_sealed = ReservoirProperties(0.3, 0.2, 0.1, Inf, 5.0)
# very small threshold (~0.18m)
const rp_quick = ReservoirProperties(0.3, 0.2, 0.1, 1000.0, 5.0)

# ---------------------------------------------------------------------------
# Injection events
# ---------------------------------------------------------------------------

# Inject at a sigle location
injection_events = Vector{Vector{InjectionEvent}}()
for i in 1:2
    pad = boundary_condition == :closed ? pad_width : 0
    topo_size = (TEST_NX + 2 * pad, TEST_NY + 2 * pad)
    if i == 1
        rate = zeros(topo_size)
        rate[div(TEST_NX, 2)+pad, div(TEST_NY, 2)+pad] = TOTAL_INJECTION_RATE
        push!(injection_events, [
            InjectionEvent(INJECTION_START_T, rate),
            InjectionEvent(INJECTION_END_T, zeros(topo_size)),
        ])
    else
        push!(injection_events, [InjectionEvent(INJECTION_START_T, zeros(topo_size))])
    end
end

# Inject at two locations at the same time.
dual_injection_events = Vector{Vector{InjectionEvent}}()
for i in 1:2
    pad = boundary_condition == :closed ? pad_width : 0
    topo_size = (TEST_NX + 2 * pad, TEST_NY + 2 * pad)
    if i == 1
        rate = zeros(topo_size)
        rate[1+pad, 1+pad] = TOTAL_INJECTION_RATE / 2
        rate[TEST_NX+pad, TEST_NY+pad] = TOTAL_INJECTION_RATE / 2
        push!(dual_injection_events, [
            InjectionEvent(INJECTION_START_T, rate),
            InjectionEvent(INJECTION_END_T, zeros(topo_size)),
        ])
    else
        push!(dual_injection_events, [InjectionEvent(INJECTION_START_T, zeros(topo_size))])
    end
end
injection_events = dual_injection_events


@testset "CO2BatchFill.jl" begin
    include("test_structs.jl")
    include("test_unit_conversion.jl")
    include("test_layer_analysis.jl")
    include("test_utils.jl")
    include("test_fill_layers.jl")
    include("test_leakage.jl")
    include("test_fill_layer.jl")
    include("test_analysis.jl")
end
