import SurfaceWaterIntegratedModeling as SWIM

@testset "Leakage unit tests" begin

    tstruct    = DOME_SCENARIO.layers[1].trap_structure
    z_vol_tables = SWIM._compute_z_vol_tables(tstruct)
    n_traps    = SWIM.numtraps(tstruct)

    @testset "find_leakage_location" begin
        # Should return a valid CartesianIndex within the padded topography bounds
        loc = find_leakage_location(1, tstruct)
        @test loc isa CartesianIndex{2}
        nx, ny = size(tstruct.topography)
        @test 1 <= loc[1] <= nx
        @test 1 <= loc[2] <= ny

        # The leakage location should be at the minimum depth in the trap's footprint
        footprint  = tstruct.footprints[1]
        topo_vals  = tstruct.topography[footprint]
        @test tstruct.topography[loc] == minimum(topo_vals)
    end

    @testset "compute_leakage_volume - sealed caprock" begin
        # Infinite leakage height: trap always spills before it can leak → nothing
        for trap_id in 1:n_traps
            result = compute_leakage_volume(trap_id, z_vol_tables[trap_id], tstruct, Inf)
            @test isnothing(result)
        end
    end

    @testset "compute_leakage_volume - leaky caprock (rp_quick)" begin
        # Very small leakage height (~0.18m): leakage elevation is well below spillpoint,
        # so leakage volume should be a small but non-negative finite value.
        h = rp_quick.leakage_height
        @test isfinite(h)
        @test h > 0.0

        result = compute_leakage_volume(1, z_vol_tables[1], tstruct, h)
        # Leakage elevation is far below the spillpoint, so result should not be nothing
        @test !isnothing(result)
        @test result >= 0.0
        # Leakage volume must be less than the total trap capacity
        @test result < maximum(z_vol_tables[1][2])
    end

    @testset "Sealed layer: initialize_leakage_state" begin
        state = initialize_leakage_state(
            tstruct, z_vol_tables, rp_sealed,
            rp_sealed.sand_residual_co2_saturation,
            rp_sealed.residual_leakage_time
        )
        @test length(state.leaking)  == n_traps
        @test length(state.draining) == n_traps
        @test length(state.leakage_height) == n_traps

        @test all(.!state.leaking)
        @test all(.!state.draining)
        @test all(state.leakage_start_time .== Inf)
        @test all(state.leakage_height .== Inf)
    end

    @testset "Leaky layer: initialize_leakage_state" begin
        state = initialize_leakage_state(
            tstruct, z_vol_tables, rp_quick,
            rp_quick.sand_residual_co2_saturation,
            rp_quick.residual_leakage_time
        )
        @test length(state.leaking) == n_traps
        @test all(.!state.leaking)  # nothing leaking at t=0

        @test all(isfinite.(state.leakage_height))
        @test all(state.leakage_height .< 1.0)
        @test all(state.leakage_volume .>= 0.0)
        @test all(isfinite.(state.leakage_volume))
    end

    @testset "GRF topology: initialize_leakage_state" begin
        tstruct_grf  = GRF_SCENARIO.layers[1].trap_structure
        z_vol_grf    = SWIM._compute_z_vol_tables(tstruct_grf)
        n_grf        = SWIM.numtraps(tstruct_grf)

        state = initialize_leakage_state(
            tstruct_grf, z_vol_grf, rp_quick,
            rp_quick.sand_residual_co2_saturation,
            rp_quick.residual_leakage_time
        )
        @test length(state.leaking)         == n_grf
        @test length(state.leakage_volume)  == n_grf
        @test all(.!state.leaking)
        @test all(state.leakage_start_time .== Inf)
    end

end
