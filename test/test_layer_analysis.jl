# Tests for src/layer_analysis.jl

@testset "Topography interface and layer analysis" begin
    topo_local = dome_topo

    @test get_grid_dimensions(topo_local) == (TEST_NX, TEST_NY)
    @test get_grid_spacing(topo_local) == (TEST_DX, TEST_DY)
    @test get_num_layers(topo_local) == 2

    domain_local = create_domain(topo_local, 1.0)
    @test domain_local.nx == TEST_NX
    @test domain_local.ny == TEST_NY

    # Analyze layers with closed boundary conditions
    layers_local = analyze_base_surfaces(topo_local; boundary_condition=:closed)
    @test length(layers_local) == 2
    @test layers_local[1].boundary_condition == :closed
    for layer in layers_local
        @test SurfaceWaterIntegratedModeling.numtraps(layer.trap_structure) > 0
    end

    # Analyze layers with open boundary conditions
    layers_open = analyze_base_surfaces(topo_local; boundary_condition=:open)
    @test length(layers_open) == 2
    @test layers_open[1].boundary_condition == :open
    for layer in layers_open
        @test SurfaceWaterIntegratedModeling.numtraps(layer.trap_structure) > 0
    end
end

@testset "GRF topology" begin
    n_grf_traps = SurfaceWaterIntegratedModeling.numtraps(grf_layers[1].trap_structure)
    @test n_grf_traps > 1

    # Both domes have the same grid size, so dimensions are preserved
    @test get_grid_dimensions(grf_topo) == (TEST_NX, TEST_NY)
    @test get_num_layers(grf_topo) == 2
end
