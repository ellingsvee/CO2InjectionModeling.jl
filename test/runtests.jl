using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using Random
using Test

Random.seed!(1)

# Test grid parameters
const TEST_NX, TEST_NY = 100, 100
const TEST_LENGTH_X, TEST_LENGTH_Y = 1000.0, 1000.0
const TEST_DX, TEST_DY = TEST_LENGTH_X / TEST_NX, TEST_LENGTH_Y / TEST_NY
const _thick = 10.0
const boundary_condition = :closed  # :open or :closed
const pad_width = 2
const RESIDUAL_TRAPPING = 0.4

const TOTAL_INJECTION_RATE = 20_000.0
const INJECTION_START_T = 0.0
const INJECTION_END_T = 15.0

const MASS_ATOL = 1e-8

# Scenario abstraction
struct TestScenario
    name::String
    topo::GenericTopography
    domain::Domain3D
    layers::Vector              # Vector{Layer} from analyze_base_surfaces
    injection_events::Vector{Vector{InjectionEvent}}   # per-layer
end

# Scenario factory helpers
function _make_dome(nx, ny, dx, dy, base_depth;
    dome_center=(0.5, 0.5), dome_sigma=150.0, dome_amp=4.0)
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
Build dual-corner injection events for a two-layer system.
Layer 1 (deepest): inject at two opposite corners; Layer 2: zero injection.
"""
function _make_dual_injection_events(; nx=TEST_NX, ny=TEST_NY,
    bc=boundary_condition, pw=pad_width,
    total_rate=TOTAL_INJECTION_RATE,
    t_start=INJECTION_START_T, t_end=INJECTION_END_T)
    events = Vector{Vector{InjectionEvent}}()
    p = bc == :closed ? pw : 0
    topo_size = (nx + 2 * p, ny + 2 * p)
    for i in 1:2
        if i == 1
            rate = zeros(topo_size)
            rate[1+p, 1+p] = total_rate / 2
            rate[nx+p, ny+p] = total_rate / 2
            push!(events, [
                InjectionEvent(t_start, rate),
                InjectionEvent(t_end, zeros(topo_size)),
            ])
        else
            push!(events, [InjectionEvent(t_start, zeros(topo_size))])
        end
    end
    return events
end

"""
Build a TestScenario from two raw elevation surfaces (top = shallowest, bottom = deepest).
"""
function _assemble_scenario(name::String, top_surface::Matrix, bottom_surface::Matrix, injection_events::Vector{Vector{InjectionEvent}})
    # sand_layers ordered shallowest-first, deepest-last
    _layers = [
        Dict{String,Any}("name" => "L2", "top" => top_surface, "base" => top_surface .+ _thick),
        Dict{String,Any}("name" => "L1", "top" => bottom_surface, "base" => bottom_surface .+ _thick),
    ]
    topo = GenericTopography(_layers, TEST_NX, TEST_NY, TEST_DX, TEST_DY,
        minimum(top_surface), maximum(bottom_surface) + _thick)
    domain = create_domain(topo, 1.0)
    analyzed_layers = analyze_base_surfaces(topo; boundary_condition=boundary_condition, pad_width=pad_width)
    return TestScenario(name, topo, domain, analyzed_layers, injection_events)
end

function make_dome_scenario(
    injection_events::Vector{Vector{InjectionEvent}}=_make_dual_injection_events()
)
    top = _make_dome(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 850.0)
    bottom = _make_dome(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 900.0)
    return _assemble_scenario("dome", top, bottom, injection_events)
end

function make_grf_scenario(
    injection_events::Vector{Vector{InjectionEvent}}=_make_dual_injection_events()
)
    cov = CovarianceFunction(2, Matern(200, 2, σ=3.0))
    pts_x = range(0, stop=(TEST_NX - 1) * TEST_DX, step=TEST_DX)
    pts_y = range(0, stop=(TEST_NY - 1) * TEST_DY, step=TEST_DY)
    top = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y, minpadding=113)) .+ 900.0
    bottom = sample(GaussianRandomField(cov, CirculantEmbedding(), pts_x, pts_y, minpadding=113)) .+ 850.0
    return _assemble_scenario("grf", top, bottom, injection_events)
end

# Build scenarios (at module level so included test files can access them)
const DOME_SCENARIO = make_dome_scenario()
const GRF_SCENARIO = make_grf_scenario()
const ALL_SCENARIOS = [DOME_SCENARIO, GRF_SCENARIO]

# Reservoir properties (shared across all tests)
const rp_sealed = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, Inf, 5.0)
# very small leakage threshold (~0.18 m) — forces leakage in reasonable simulation time
const rp_quick = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, 1000.0, 5.0)

function rp_non_stat(tstruct::TrapStructure)
    # Initialise a vector of random capillary entry pressures for each trap
    num_traps = numtraps(tstruct)
    random_pressures = rand(num_traps) .* 5000.0 .+ 500.0  # Random pressures between 500 Pa and 5500 Pa
    return ReservoirProperties(
        0.3, RESIDUAL_TRAPPING, 0.1,
        random_pressures,
        5.0
    )
end




# Shared test helpers
"""
Assert that a LayerSnapshot satisfies exact mass conservation.
"""
function assert_mass_conservation(snap::LayerSnapshot; atol=MASS_ATOL)
    @test snap.total_stored >= -atol
    @test snap.total_drained >= -atol
    @test total_to_next_layer(snap) >= -atol
    @test total_passthrough(snap) >= -atol
    # stored + upward leakage == total injected
    @test snap.total_stored + total_to_next_layer(snap) ≈ snap.total_injected atol = atol
end

"""
Compute expected total CO2 injected (in SWIM units) up to time `t` directly from
the InjectionEvent schedule, without going through the WeatherEvent machinery.
This serves as an independent ground-truth reference for `snap.total_injected`.
"""
function total_injected_from_schedule(
    injection_events::Vector{InjectionEvent},
    rp::ReservoirProperties, domain::Domain3D, t::Float64
)::Float64
    total_physical = 0.0
    n = length(injection_events)
    for i in 1:n
        t_start = injection_events[i].timestamp
        rate = sum(injection_events[i].injection_rate)
        t_end = i < n ? min(t, injection_events[i+1].timestamp) : t
        total_physical += rate * max(0.0, t_end - t_start)
    end
    return physical_volume_to_swim_volume(total_physical, rp, domain)
end

# Run all tests
@testset "CO2BatchFill.jl" begin
    # Non-parameterized tests (scenario-independent)
    include("test_structs.jl")
    include("test_unit_conversion.jl")
    include("test_layer_analysis.jl")
    include("test_utils.jl")
    include("test_leakage.jl")

    # Load parameterized test functions (define run_* functions)
    include("test_fill_layer.jl")
    include("test_analysis.jl")
    include("test_fill_layers.jl")
    include("test_multilayer_analysis.jl")

    # Run parameterized tests against every scenario
    for scenario in ALL_SCENARIOS
        @testset "Scenario: $(scenario.name)" begin
            run_fill_layer_tests(scenario)
            run_analysis_tests(scenario)
            run_fill_layers_tests(scenario)
            run_multilayer_analysis_tests(scenario)
        end
    end
end
