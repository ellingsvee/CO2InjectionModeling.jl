using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using Random
using Test

Random.seed!(1)

# # ---------------------------------------------------------------------------
# # Test grid parameters 
# # ---------------------------------------------------------------------------
const TEST_NX, TEST_NY = 20, 20
const TEST_LENGTH_X, TEST_LENGTH_Y = 1000.0, 1000.0
const TEST_DX, TEST_DY = TEST_LENGTH_X / TEST_NX, TEST_LENGTH_Y / TEST_NY
const _thick = 10.0
const boundary_condition = :closed # :open or :closed


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
    analyzed_layers = analyze_base_surfaces(topo; boundary_condition=boundary_condition)
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
pts_x = range(0, stop=TEST_NX * (TEST_DX - 1), step=TEST_DX)
pts_y = range(0, stop=TEST_NY * (TEST_DY - 1), step=TEST_DY)

const _grf1 = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y)) .+ 900.0
const _grf2 = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y)) .+ 850.0

grf_topo, grf_domain, grf_layers = create_domain_and_layers(_grf1, _grf2, _thick, TEST_NX, TEST_NY, TEST_DX, TEST_DY, boundary_condition)

# ---------------------------------------------------------------------------
# Reservoir properties
# ---------------------------------------------------------------------------

const rp_leaky = ReservoirProperties(0.3, 0.2, 0.1, 12_000.0, 5.0)
const rp_sealed = ReservoirProperties(0.3, 0.2, 0.1, Inf, 5.0)
# rp_quick: very small threshold (~0.18m) to guarantee leakage in tests without large injection
const rp_quick = ReservoirProperties(0.3, 0.2, 0.1, 1000.0, 5.0)

# ---------------------------------------------------------------------------
# Injection events
# ---------------------------------------------------------------------------

injection_events = Vector{Vector{InjectionEvent}}()
for i in 1:2
    pad = boundary_condition == :closed ? 1 : 0
    topo_size = (TEST_NX + 2 * pad, TEST_NY + 2 * pad)
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
