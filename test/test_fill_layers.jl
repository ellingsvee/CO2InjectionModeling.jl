import SurfaceWaterIntegratedModeling as SWIM

function run_fill_layers_tests(scenario::TestScenario)
    @testset "Multi-layer leakage weather events" begin

        layer_idx = 1
        layer1    = scenario.layers[layer_idx]
        we1       = convert_injection_event_to_weather_event(
            scenario.injection_events[layer_idx], rp_quick, scenario.domain)
        seq1, ls1 = fill_sequence_with_leakage(
            layer1.trap_structure, rp_quick, we1; verbose=false)

        we2_direct = convert_injection_event_to_weather_event(
            scenario.injection_events[2], rp_quick, scenario.domain)


        @testset "No leakage -> returns direct events unchanged" begin
            we1_sealed   = convert_injection_event_to_weather_event(
                scenario.injection_events[layer_idx], rp_sealed, scenario.domain)
            seq_sealed, ls_sealed = fill_sequence_with_leakage(
                layer1.trap_structure, rp_sealed, we1_sealed; verbose=false)

            result = generate_leakage_weather_events(
                seq_sealed, ls_sealed, layer1.trap_structure,
                rp_sealed, rp_sealed, we2_direct)

            @test length(result) == length(we2_direct)
            for (r, d) in zip(result, we2_direct)
                @test r.timestamp == d.timestamp
                if r.rain_rate isa Matrix && d.rain_rate isa Matrix
                    @test r.rain_rate ≈ d.rain_rate atol = MASS_ATOL
                else
                    @test r.rain_rate == d.rain_rate
                end
            end
        end


        @testset "Leaking layer -> non-zero events generated" begin
            @test length(ls1.leakage_records) > 0

            we2_combined = generate_leakage_weather_events(
                seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

            @test !isempty(we2_combined)
            @test any(any(we.rain_rate .> 0) for we in we2_combined)

            ts = [we.timestamp for we in we2_combined]
            @test issorted(ts)
            # All weather-event rain rates must be non-negative (no ghost CO2)
            @test all(all(we.rain_rate .>= 0) for we in we2_combined)

            @test we2_combined[1].timestamp <= seq1[1].timestamp + MASS_ATOL
        end


        @testset "Mass conservation: layer-1 output = layer-2 input" begin
            @test length(ls1.leakage_records) > 0

            we2_combined = generate_leakage_weather_events(
                seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

            t_sim_start  = seq1[1].timestamp
            t_leak_start = ls1.leakage_start_time[findfirst(ls1.leaking)]
            t_drain_end  = t_leak_start + rp_quick.residual_leakage_time

            # 8 timepoints spanning before-leak → during-drain → after-drain
            timepoints = sort(unique(filter(isfinite, [
                t_sim_start,
                t_leak_start,
                t_leak_start + 0.25 * rp_quick.residual_leakage_time,
                t_leak_start + 0.5  * rp_quick.residual_leakage_time,
                t_drain_end,
                t_drain_end + 0.5 * rp_quick.residual_leakage_time,
                t_drain_end + rp_quick.residual_leakage_time,
                t_drain_end + 2.0 * rp_quick.residual_leakage_time,
            ])))

            snaps1 = generate_layer_snapshots(layer1, layer_idx, seq1, ls1, we1, timepoints)

            for (snap, t) in zip(snaps1, timepoints)
                left_layer1    = total_upward_leakage(snap)
                received_layer2 = compute_total_injected(we2_combined, t) -
                                  compute_total_injected(we2_direct, t)
                @test left_layer1 ≈ received_layer2 atol = MASS_ATOL
            end

            # total_to_next_layer is non-decreasing
            for i in 2:length(snaps1)
                @test total_to_next_layer(snaps1[i]) >= total_to_next_layer(snaps1[i-1]) - MASS_ATOL
            end
        end


        @testset "Drainage component matches residual drainage" begin
            we2_combined = generate_leakage_weather_events(
                seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

            t_leak_start = ls1.leakage_start_time[findfirst(ls1.leaking)]
            t_drain_end  = t_leak_start + rp_quick.residual_leakage_time
            t_after      = t_drain_end + 1.0

            total_drained_layer1 = compute_total_drained(t_after, ls1)
            direct_total         = compute_total_injected(we2_direct,   t_after)
            leakage_total        = compute_total_injected(we2_combined, t_after) - direct_total

            @test leakage_total >= total_drained_layer1 - MASS_ATOL

            snaps1 = generate_layer_snapshots(layer1, layer_idx, seq1, ls1, we1, [t_after])
            @test total_upward_leakage(snaps1[1]) ≈ leakage_total atol = MASS_ATOL
        end


        @testset "End-to-end: layer-2 simulation runs with generated events" begin
            we2_combined = generate_leakage_weather_events(
                seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

            layer2 = scenario.layers[2]   # upper layer (fixed: was dome_layers[2])

            seq2, ls2 = fill_sequence_with_leakage(
                layer2.trap_structure, rp_quick, we2_combined; verbose=false)

            @test length(seq2) > 0
            @test all(isfinite(se.timestamp) for se in seq2)

            # Mass conservation within layer 2 at a late timepoint
            t_check = seq2[end].timestamp
            snap2   = generate_layer_snapshot(
                layer2, 2, seq2, ls2, we2_combined, t_check,
                SWIM.trap_states_at_timepoints(layer2.trap_structure, seq2, [t_check];
                                               verbose=false)[1])

            @test snap2.total_stored         >= -MASS_ATOL
            @test snap2.total_drained        >= -MASS_ATOL
            @test total_to_next_layer(snap2) >= -MASS_ATOL
            @test snap2.total_stored + total_to_next_layer(snap2) ≈ snap2.total_injected atol = MASS_ATOL
        end

    end
end
