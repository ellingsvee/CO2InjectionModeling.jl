using SurfaceWaterIntegratedModeling: numtraps

@testset "Fill sequence for layer" begin

    @testset "Sealed caprock - single dome" begin
        # Sanity check: fill runs and returns a non-empty sequence
        layer_idx = 1
        tstruct = dome_layers[layer_idx].trap_structure
        we_events = convert_injection_event_to_weather_event(injection_events[layer_idx], rp_sealed, dome_domain)
        seq, _ = fill_sequence_with_leakage(tstruct, rp_sealed, we_events; verbose=false)
        @test length(seq) > 0
    end


    @testset "Two-dome topology - multi-trap fill" begin
        # Two distinct structural highs should create at least 2 traps
        layer_idx = 1
        tstruct = grf_layers[layer_idx].trap_structure
        @test numtraps(tstruct) >= 2

        # Fill sequence should run without error on multi-trap topology
        we_events = convert_injection_event_to_weather_event(injection_events[layer_idx], rp_sealed, grf_domain)
        seq, _ = fill_sequence_with_leakage(tstruct, rp_sealed, we_events)
        @test length(seq) > 0

        # With injection at the left dome and sealed caprock, CO2 fills traps in spillpoint order
        # (at minimum the injected trap should fill and spill)
        fill_times = [e.timestamp for e in seq]
        @test all(isfinite, fill_times)
    end


    @testset "Leaking layer" begin
        # Sanity check: fill runs and returns a non-empty sequence
        layer_idx = 1
        tstruct = dome_layers[layer_idx].trap_structure
        we_events = convert_injection_event_to_weather_event(injection_events[layer_idx], rp_quick, dome_domain)
        seq, leakage_state = fill_sequence_with_leakage(tstruct, rp_quick, we_events; verbose=false)

        @test length(seq) > 0
        @test length(leakage_state.leakage_records) > 0

    end

end
