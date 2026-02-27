using CO2BatchFill
using SurfaceWaterIntegratedModeling
using Test

# ---------------------------------------------------------------------------
# Helper: create a depth surface with regional dip + one Gaussian dome
# ---------------------------------------------------------------------------
function make_test_surface(nx, ny, dx, dy, base_depth;
                           dome_center=(0.5, 0.5),  # fraction of domain
                           dome_sigma=150.0, dome_amp=4.0,
                           slope_x=0.002, slope_y=0.001)
    xs = [(i - 0.5) * dx for i in 1:nx]
    ys = [(j - 0.5) * dy for j in 1:ny]
    cx = dome_center[1] * nx * dx
    cy = dome_center[2] * ny * dy
    surface = zeros(nx, ny)
    for i in 1:nx, j in 1:ny
        surface[i, j] = base_depth + slope_x * xs[i] + slope_y * ys[j]
        r2 = (xs[i] - cx)^2 + (ys[j] - cy)^2
        surface[i, j] -= dome_amp * exp(-r2 / (2 * dome_sigma^2))
    end
    return surface
end

# ---------------------------------------------------------------------------
# Test grid parameters (kept small for speed)
# ---------------------------------------------------------------------------
const TEST_NX, TEST_NY = 20, 20
const TEST_DX, TEST_DY = 50.0, 50.0

@testset "CO2BatchFill.jl" begin

    # -----------------------------------------------------------------------
    @testset "Data structures" begin
        domain = Domain3D(10, 10, 5, 100.0, 100.0, 50.0, 800.0, 850.0)
        @test domain.nx == 10
        @test domain.dx ≈ 10.0
        @test domain.dy ≈ 10.0
        @test domain.dz ≈ 10.0

        rp = ReservoirProperties(0.3, 0.2, 0.1, 50000.0, 5.0)
        @test rp.sand_porosity == 0.3
        @test isfinite(rp.leakage_height)
        @test rp.leakage_height > 0

        rp_sealed = ReservoirProperties(0.3, 0.2, 0.1, Inf, 5.0)
        @test rp_sealed.leakage_height == Inf
    end

    # -----------------------------------------------------------------------
    @testset "Volume conversion round-trip" begin
        domain = Domain3D(10, 10, 5, 100.0, 100.0, 50.0, 800.0, 850.0)
        rp = ReservoirProperties(0.3, 0.2, 0.1, 50000.0, 5.0)

        vol_physical = 1000.0
        vol_swim = physical_volume_to_swim_volume(vol_physical, rp, domain)
        vol_back = swim_volume_to_physical_volume(vol_swim, rp, domain)
        @test vol_back ≈ vol_physical
    end

    # -----------------------------------------------------------------------
    @testset "Topography interface and layer analysis" begin
        s1 = make_test_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 900.0)
        s2 = make_test_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 850.0;
                               dome_center=(0.4, 0.4))
        thick = 10.0

        # Shallowest first, deepest last (analyze_base_surfaces reverses)
        sand_layers = [
            Dict{String,Any}("name" => "L2", "top" => s2, "base" => s2 .+ thick),
            Dict{String,Any}("name" => "L1", "top" => s1, "base" => s1 .+ thick),
        ]

        topo = GenericTopography(sand_layers, TEST_NX, TEST_NY, TEST_DX, TEST_DY,
                                 minimum(s2), maximum(s1 .+ thick))

        @test get_grid_dimensions(topo) == (TEST_NX, TEST_NY)
        @test get_grid_spacing(topo) == (TEST_DX, TEST_DY)
        @test get_num_layers(topo) == 2

        domain = create_domain(topo, 1.0)
        @test domain.nx == TEST_NX
        @test domain.ny == TEST_NY

        layers = analyze_base_surfaces(topo; boundary_condition=:closed)
        @test length(layers) == 2
        @test layers[1].boundary_padding == 1
        for layer in layers
            @test SurfaceWaterIntegratedModeling.numtraps(layer.trap_structure) > 0
        end
    end

    # -----------------------------------------------------------------------
    # Shared setup for simulation tests
    # -----------------------------------------------------------------------
    s1 = make_test_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 900.0)
    s2 = make_test_surface(TEST_NX, TEST_NY, TEST_DX, TEST_DY, 850.0;
                           dome_center=(0.4, 0.4))
    thick = 10.0

    sand_layers = [
        Dict{String,Any}("name" => "L2", "top" => s2, "base" => s2 .+ thick),
        Dict{String,Any}("name" => "L1", "top" => s1, "base" => s1 .+ thick),
    ]
    topo = GenericTopography(sand_layers, TEST_NX, TEST_NY, TEST_DX, TEST_DY,
                             minimum(s2), maximum(s1 .+ thick))

    domain = create_domain(topo, 1.0)
    layers = analyze_base_surfaces(topo; boundary_condition=:closed)

    rp_leaky = ReservoirProperties(0.3, 0.2, 0.1, 50000.0, 5.0)
    rp_sealed = ReservoirProperties(0.3, 0.2, 0.1, Inf, 5.0)

    # -----------------------------------------------------------------------
    @testset "Single-layer simulation" begin
        topo_size = size(layers[1].trap_structure.topography)
        pad = layers[1].boundary_padding

        injection_rate = zeros(topo_size)
        injection_rate[div(TEST_NX, 2) + pad, div(TEST_NY, 2) + pad] = 2000.0

        weather_events = convert_injection_event_to_weather_event(
            [InjectionEvent(0.0, injection_rate)], rp_leaky, domain
        )

        tstruct = layers[1].trap_structure
        seq, leakage_state = fill_layer(tstruct, domain, rp_leaky, weather_events)

        @test !isempty(seq)
        @test length(leakage_state.leaking) == numtraps(tstruct)
    end

    # -----------------------------------------------------------------------
    @testset "Multi-layer simulation" begin
        injection_events = Vector{Vector{InjectionEvent}}()
        for (i, layer) in enumerate(layers)
            topo_size = size(layer.trap_structure.topography)
            pad = layer.boundary_padding
            if i == 1
                rate = zeros(topo_size)
                rate[div(TEST_NX, 2) + pad, div(TEST_NY, 2) + pad] = 2000.0
                push!(injection_events, [
                    InjectionEvent(0.0, rate),
                    InjectionEvent(8.0, zeros(topo_size)),
                ])
            else
                push!(injection_events, [InjectionEvent(0.0, zeros(topo_size))])
            end
        end

        rprops = [rp_leaky, rp_sealed]
        seqs, leakage_states = fill_layers(layers, domain, rprops, injection_events)

        @test length(seqs) == 2
        @test length(leakage_states) == 2
        @test !isempty(seqs[1])
    end

    # -----------------------------------------------------------------------
    @testset "Snapshot generation and mass balance" begin
        injection_events = Vector{Vector{InjectionEvent}}()
        for (i, layer) in enumerate(layers)
            topo_size = size(layer.trap_structure.topography)
            pad = layer.boundary_padding
            if i == 1
                rate = zeros(topo_size)
                rate[div(TEST_NX, 2) + pad, div(TEST_NY, 2) + pad] = 2000.0
                push!(injection_events, [
                    InjectionEvent(0.0, rate),
                    InjectionEvent(8.0, zeros(topo_size)),
                ])
            else
                push!(injection_events, [InjectionEvent(0.0, zeros(topo_size))])
            end
        end

        rprops = [rp_leaky, rp_sealed]
        seqs, leakage_states = fill_layers(layers, domain, rprops, injection_events)

        snapshots = generate_reservoir_snapshots(
            layers, seqs, leakage_states, domain, rprops, injection_events;
            num_snapshots=5, start_time=0.0, end_time=10.0, verbose=false
        )

        @test length(snapshots) == 5
        @test snapshots[1].timestamp ≈ 0.0
        @test snapshots[end].timestamp ≈ 10.0

        # After injection, there should be CO2 in the reservoir
        final = snapshots[end]
        @test final.total_injected_m3 > 0
        @test final.total_stored_m3 ≥ 0

        # Mass balance: injected = stored + leaked (within tolerance)
        @test final.mass_balance_error_percent < 1.0
    end
end
