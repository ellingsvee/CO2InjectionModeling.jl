using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using Random
using Test

Random.seed!(1)

# ---------------------------------------------------------------------------
# Helper: create a depth surface with regional dip + one Gaussian dome
# ---------------------------------------------------------------------------
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
