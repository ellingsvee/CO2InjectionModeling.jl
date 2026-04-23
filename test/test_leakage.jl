import SurfaceWaterIntegratedModeling as SWIM

@testset "Leakage unit tests" begin

    tstruct = DOME_SCENARIO.layers[1].trap_structure
    z_vol_tables = SWIM._compute_z_vol_tables(tstruct)
    n_traps = SWIM.numtraps(tstruct)

    @testset "find_leakage_location" begin
        # Should return a valid CartesianIndex within the padded topography bounds
        loc = find_leakage_location(1, tstruct)
        @test loc isa CartesianIndex{2}
        nx, ny = size(tstruct.topography)
        @test 1 <= loc[1] <= nx
        @test 1 <= loc[2] <= ny

        # The leakage location should be at the minimum depth in the trap's footprint
        footprint = tstruct.footprints[1]
        topo_vals = tstruct.topography[footprint]
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
        @test length(state.leaking) == n_traps
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

    @testset "Per-trap leakage heights (rp_non_stat)" begin
        rp_nt = rp_non_stat(tstruct)

        @test rp_nt.leakage_height isa Vector{Float64}
        @test length(rp_nt.leakage_height) == n_traps
        @test all(isfinite.(rp_nt.leakage_height))
        @test all(rp_nt.leakage_height .> 0.0)

        state = initialize_leakage_state(
            tstruct, z_vol_tables, rp_nt,
            rp_nt.sand_residual_co2_saturation,
            rp_nt.residual_leakage_time
        )

        @test length(state.leakage_height) == n_traps
        @test state.leakage_height == rp_nt.leakage_height
        @test all(.!state.leaking)

        # Each trap's leakage volume must match its own height
        for trap_id in 1:n_traps
            expected_vol = compute_leakage_volume(
                trap_id, z_vol_tables[trap_id], tstruct, rp_nt.leakage_height[trap_id]
            )
            if isnothing(expected_vol)
                @test state.leakage_volume[trap_id] == Inf
            else
                @test state.leakage_volume[trap_id] ≈ expected_vol
            end
        end

        # With multiple traps, per-trap heights should differ (random pressures)
        if n_traps > 1
            @test length(unique(state.leakage_height)) > 1
        end
    end

    @testset "Per-trap leakage heights: mixed finite and Inf" begin
        pressures = fill(1000.0, n_traps)
        for i in 2:2:n_traps
            pressures[i] = Inf
        end
        rp_mixed = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, pressures, 5.0)

        @test rp_mixed.leakage_height isa Vector{Float64}
        for i in 1:n_traps
            if i % 2 == 0
                @test rp_mixed.leakage_height[i] == Inf
            else
                @test isfinite(rp_mixed.leakage_height[i])
            end
        end

        state = initialize_leakage_state(
            tstruct, z_vol_tables, rp_mixed,
            rp_mixed.sand_residual_co2_saturation,
            rp_mixed.residual_leakage_time
        )

        for i in 2:2:n_traps
            @test state.leakage_volume[i] == Inf
        end
        for i in 1:2:n_traps
            @test isfinite(state.leakage_volume[i])
            @test state.leakage_volume[i] >= 0.0
        end
    end

    @testset "Shale-break child does not reduce parent leakage volume" begin
        # When a child trap is a shale break (pressure=0, leakage_height=0),
        # the parent's leakage volume should be computed from the effective
        # bottom excluding the break child's footprint, not from the true
        # bottom that includes the break child's deep area.
        for scenario in ALL_SCENARIOS
            tstruct_s = scenario.layers[1].trap_structure
            z_vol_s = SWIM._compute_z_vol_tables(tstruct_s)
            n = SWIM.numtraps(tstruct_s)

            # Find a parent-child pair
            parent_trap = nothing
            break_child = nothing
            for t in 1:n
                children = SWIM.subtrapsof(tstruct_s, t)
                if !isempty(children)
                    parent_trap = t
                    break_child = children[1]
                    break
                end
            end
            isnothing(parent_trap) && continue

            # Set only the child to shale-break; parent gets a normal threshold
            pressures = fill(1000.0, n)
            pressures[break_child] = 0.0
            rp_break = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, pressures, 5.0)

            # Compute leakage state with and without the break child
            state_with_break = initialize_leakage_state(
                tstruct_s, z_vol_s, rp_break,
                rp_break.sand_residual_co2_saturation,
                rp_break.residual_leakage_time
            )

            # Reference: same parent threshold but NO break children
            pressures_nobreak = fill(1000.0, n)
            rp_nobreak = ReservoirProperties(0.3, RESIDUAL_TRAPPING, 0.1, pressures_nobreak, 5.0)
            state_no_break = initialize_leakage_state(
                tstruct_s, z_vol_s, rp_nobreak,
                rp_nobreak.sand_residual_co2_saturation,
                rp_nobreak.residual_leakage_time
            )

            # Break child must have leakage_volume=0
            @test state_with_break.leakage_volume[break_child] == 0.0

            # Parent's leakage volume with a break child should be >= the
            # volume without a break child (the corrected bottom is higher,
            # so leakage_elevation is higher, requiring more volume)
            parent_capacity = tstruct_s.trapvolumes[parent_trap] - tstruct_s.subvolumes[parent_trap]
            if parent_capacity > 0
                @test state_with_break.leakage_volume[parent_trap] >= state_no_break.leakage_volume[parent_trap]
            end
        end
    end

    @testset "GRF topology: initialize_leakage_state" begin
        tstruct_grf = GRF_SCENARIO.layers[1].trap_structure
        z_vol_grf = SWIM._compute_z_vol_tables(tstruct_grf)
        n_grf = SWIM.numtraps(tstruct_grf)

        state = initialize_leakage_state(
            tstruct_grf, z_vol_grf, rp_quick,
            rp_quick.sand_residual_co2_saturation,
            rp_quick.residual_leakage_time
        )
        @test length(state.leaking) == n_grf
        @test length(state.leakage_volume) == n_grf
        @test all(.!state.leaking)
        @test all(state.leakage_start_time .== Inf)
    end

end
