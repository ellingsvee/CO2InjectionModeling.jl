using SurfaceWaterIntegratedModeling: numtraps

function run_analysis_tests(scenario::TestScenario)
    @testset "Layer snapshots and mass conservation" begin

        @testset "Sealed caprock – no leakage" begin
            layer_idx = 1
            layer     = scenario.layers[layer_idx]
            inj       = scenario.injection_events[layer_idx]
            we_events = convert_injection_event_to_weather_event(inj, rp_sealed, scenario.domain)

            seq, leakage_state = fill_sequence_with_leakage(
                layer.trap_structure, rp_sealed, we_events; verbose=false)

            t_sim_start = seq[1].timestamp
            # Four timepoints: start, mid-injection, end-injection, 2× injection duration
            t_inj_dur   = INJECTION_END_T - INJECTION_START_T
            timepoints  = sort(unique([
                t_sim_start,
                t_sim_start + 0.5 * t_inj_dur,
                t_sim_start + t_inj_dur,
                t_sim_start + 2.0 * t_inj_dur,
            ]))

            snaps = generate_layer_snapshots(layer, layer_idx, seq, leakage_state,
                we_events, timepoints)

            for snap in snaps
                assert_mass_conservation(snap)
            end

            # With sealed caprock nothing drains or leaves the layer
            for snap in snaps
                @test snap.total_drained        ≈ 0.0 atol = MASS_ATOL
                @test total_to_next_layer(snap) ≈ 0.0 atol = MASS_ATOL
                @test snap.total_stored         ≈ snap.total_injected atol = MASS_ATOL
            end

            # Monotonically accumulating CO2
            @test snaps[1].total_injected <= snaps[end].total_injected + MASS_ATOL

            # Independent ground-truth check: WeatherEvent total matches injection schedule
            for snap in snaps
                expected = total_injected_from_schedule(inj, rp_sealed, scenario.domain, snap.timestamp)
                @test snap.total_injected ≈ expected atol = MASS_ATOL
            end
        end


        @testset "Leaking layer: Mass conservation with drainage" begin
            layer_idx = 1
            layer     = scenario.layers[layer_idx]
            inj       = scenario.injection_events[layer_idx]
            we_events = convert_injection_event_to_weather_event(inj, rp_quick, scenario.domain)

            seq, leakage_state = fill_sequence_with_leakage(
                layer.trap_structure, rp_quick, we_events; verbose=false)

            @test length(leakage_state.leakage_records) > 0  # leakage must have occurred

            t_sim_start  = seq[1].timestamp
            t_leak_start = leakage_state.leakage_start_time[findfirst(leakage_state.leaking)]
            t_drain_end  = t_leak_start + rp_quick.residual_leakage_time

            # 8 timepoints spanning before-leak → during-drain → after-drain
            timepoints = sort(unique([
                t_sim_start,
                t_leak_start - 0.1 * rp_quick.residual_leakage_time,
                t_leak_start,
                t_leak_start + 0.25 * rp_quick.residual_leakage_time,
                t_leak_start + 0.5  * rp_quick.residual_leakage_time,
                t_drain_end,
                t_drain_end + 0.5 * rp_quick.residual_leakage_time,
                t_drain_end + rp_quick.residual_leakage_time,
            ]))
            # Clamp so we don't go before the simulation start
            timepoints = filter(t -> t >= t_sim_start, timepoints)

            snaps = generate_layer_snapshots(layer, layer_idx, seq, leakage_state,
                we_events, timepoints)

            for snap in snaps
                assert_mass_conservation(snap)
            end

            # Monotonicity: CO2 doesn't flow back
            for i in 2:length(snaps)
                @test total_to_next_layer(snaps[i]) >= total_to_next_layer(snaps[i-1]) - MASS_ATOL
                @test snaps[i].total_drained         >= snaps[i-1].total_drained - MASS_ATOL
            end

            # Independent ground-truth check at each timepoint
            for snap in snaps
                expected = total_injected_from_schedule(inj, rp_quick, scenario.domain, snap.timestamp)
                @test snap.total_injected ≈ expected atol = MASS_ATOL
            end

            # After drainage window: residual volumes match expected residual saturation
            snap_after = snaps[end]
            for trap in 1:numtraps(layer.trap_structure)
                if leakage_state.draining[trap]
                    t_drain_end_trap = leakage_state.leakage_start_time[trap] + rp_quick.residual_leakage_time
                    if snap_after.timestamp >= t_drain_end_trap
                        expected_final = leakage_state.initial_volume_at_leak[trap] *
                                         rp_quick.sand_residual_co2_saturation
                        @test snap_after.trap_volumes[trap] ≈ expected_final atol = MASS_ATOL
                    end
                end
            end
        end


        @testset "Snapshot fields are self-consistent" begin
            layer_idx = 1
            layer     = scenario.layers[layer_idx]
            inj       = scenario.injection_events[layer_idx]
            we_events = convert_injection_event_to_weather_event(inj, rp_quick, scenario.domain)
            seq, leakage_state = fill_sequence_with_leakage(
                layer.trap_structure, rp_quick, we_events; verbose=false)

            t    = seq[1].timestamp + 8.0
            snap = generate_layer_snapshot(
                layer, layer_idx, seq, leakage_state, we_events, t,
                trap_states_at_timepoints(layer.trap_structure, seq, [t]; verbose=false)[1])

            @test snap.total_stored ≈ sum(snap.trap_volumes) atol = MASS_ATOL
            @test all(snap.trap_volumes .>= -MASS_ATOL)

            # leaking implies draining
            for trap in eachindex(snap.trap_leaking)
                if snap.trap_leaking[trap]
                    @test snap.trap_draining[trap]
                end
            end

            @test snap.layer_idx  == layer_idx
            @test snap.layer_name == layer.name
            @test snap.timestamp  == t
        end

    end

    @testset "Per-trap pressure: mass conservation with drainage" begin
        layer_idx = 1
        layer     = scenario.layers[layer_idx]
        tstruct   = layer.trap_structure
        inj       = scenario.injection_events[layer_idx]
        rp_nt     = rp_non_stat(tstruct)
        we_events = convert_injection_event_to_weather_event(inj, rp_nt, scenario.domain)

        seq, leakage_state = fill_sequence_with_leakage(
            tstruct, rp_nt, we_events; verbose=false)

        # Leakage heights should be per-trap
        @test leakage_state.leakage_height isa Vector{Float64}

        t_sim_start = seq[1].timestamp
        # Build timepoints spanning the simulation
        t_inj_dur = INJECTION_END_T - INJECTION_START_T
        timepoints = sort(unique(filter(isfinite, [
            t_sim_start,
            t_sim_start + 0.5 * t_inj_dur,
            t_sim_start + t_inj_dur,
            t_sim_start + 2.0 * t_inj_dur,
        ])))
        timepoints = filter(t -> t >= t_sim_start, timepoints)

        snaps = generate_layer_snapshots(layer, layer_idx, seq, leakage_state,
            we_events, timepoints)

        for snap in snaps
            assert_mass_conservation(snap)
        end

        # Independent ground-truth check
        for snap in snaps
            expected = total_injected_from_schedule(inj, rp_nt, scenario.domain, snap.timestamp)
            @test snap.total_injected ≈ expected atol = MASS_ATOL
        end

        # Monotonicity
        for i in 2:length(snaps)
            @test total_to_next_layer(snaps[i]) >= total_to_next_layer(snaps[i-1]) - MASS_ATOL
            @test snaps[i].total_drained >= snaps[i-1].total_drained - MASS_ATOL
        end
    end

end
