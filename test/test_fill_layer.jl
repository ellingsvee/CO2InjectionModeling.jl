using SurfaceWaterIntegratedModeling: numtraps

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
    end
end
