using SurfaceWaterIntegratedModeling: numtraps
import Graphs

function run_fill_layer_tests(scenario::TestScenario)
    @testset "Fill sequence for layer" begin

        @testset "Sealed caprock" begin
            layer_idx = 1
            tstruct = scenario.layers[layer_idx].trap_structure
            we_events = convert_injection_event_to_weather_event(
                scenario.injection_events[layer_idx], rp_sealed, scenario.domain)
            seq, _ = fill_sequence_with_leakage(tstruct, rp_sealed, we_events; verbose=false)
            @test length(seq) > 0
        end

        @testset "Multi-trap fill" begin
            # Any scenario should produce at least one trap
            layer_idx = 1
            tstruct = scenario.layers[layer_idx].trap_structure
            @test numtraps(tstruct) >= 1

            we_events = convert_injection_event_to_weather_event(
                scenario.injection_events[layer_idx], rp_sealed, scenario.domain)
            seq, _ = fill_sequence_with_leakage(tstruct, rp_sealed, we_events)
            @test length(seq) > 0

            fill_times = [e.timestamp for e in seq]
            @test all(isfinite, fill_times)
        end

        @testset "Leaking layer" begin
            layer_idx = 1
            tstruct = scenario.layers[layer_idx].trap_structure
            we_events = convert_injection_event_to_weather_event(
                scenario.injection_events[layer_idx], rp_quick, scenario.domain)
            seq, leakage_state = fill_sequence_with_leakage(
                tstruct, rp_quick, we_events; verbose=false)

            @test length(seq) > 0
            @test length(leakage_state.leakage_records) > 0
        end

        @testset "Shale-break with sibling cycle" begin
            # Verify that fill_sequence_with_leakage runs correctly when a
            # shale-break trap has siblings. The sibling cycle check in
            # update_spillgraph! must not redirect siblings away from the break.
            layer_idx = 1
            tstruct = scenario.layers[layer_idx].trap_structure
            n_traps = numtraps(tstruct)
            n_traps < 2 && return

            rp_sb = rp_shale_breaks(tstruct; shale_break_fraction=0.3)
            we_events = convert_injection_event_to_weather_event(
                scenario.injection_events[layer_idx], rp_sb, scenario.domain)
            seq, leakage_state = fill_sequence_with_leakage(
                tstruct, rp_sb, we_events; verbose=false)

            @test length(seq) > 0

            # Traps with Inf leakage height must never become leaking
            for t in 1:n_traps
                if rp_sb.leakage_height[t] == Inf
                    @test !leakage_state.leaking[t]
                end
            end
        end
    end
end
