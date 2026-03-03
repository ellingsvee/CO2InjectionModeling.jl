using SurfaceWaterIntegratedModeling: numtraps

@testset "Fill sequence for layer" begin

    @testset "Sealed caprock - single dome" begin
        # Sanity check: fill runs and returns a non-empty sequence
        layer_idx = 1
        tstruct = layers[layer_idx].trap_structure
        we_events = convert_injection_event_to_weather_event(injection_events[layer_idx], rp_sealed, domain)
        seq = fill_sequence_with_leakage(tstruct, rp_sealed, we_events; verbose=false)
        @test length(seq) > 0
    end

    # @testset "Sealed vs leaky produce different sequences" begin
    #     # With rp_quick (leakage_height ≈ 0.18m), leakage triggers almost immediately,
    #     # so the fill sequence should differ from the sealed case (traps behave differently).
    #     layer_idx = 1
    #     tstruct = layers[layer_idx].trap_structure
    #     we_sealed = convert_injection_event_to_weather_event(injection_events[layer_idx], rp_sealed, domain)
    #     we_quick = convert_injection_event_to_weather_event(injection_events[layer_idx], rp_quick, domain)

    #     seq_sealed = fill_sequence_with_leakage(tstruct, rp_sealed, we_sealed)
    #     seq_quick = fill_sequence_with_leakage(tstruct, rp_quick, we_quick)

    #     # Both should produce valid (non-empty) sequences
    #     @test length(seq_sealed) > 0
    #     @test length(seq_quick) > 0

    #     # The sequences should differ: leakage means CO2 drains, producing a different event timeline
    #     timestamps_sealed = [e.timestamp for e in seq_sealed]
    #     timestamps_quick = [e.timestamp for e in seq_quick]
    #     @test timestamps_sealed != timestamps_quick
    # end

    @testset "Two-dome topology - multi-trap fill" begin
        # Two distinct structural highs should create at least 2 traps
        layer_idx = 1
        tstruct = td_layers[layer_idx].trap_structure
        @test numtraps(tstruct) >= 2

        # Fill sequence should run without error on multi-trap topology
        we_events = convert_injection_event_to_weather_event(td_injection_events[layer_idx], rp_sealed, domain)
        seq = fill_sequence_with_leakage(tstruct, rp_sealed, we_events)
        @test length(seq) > 0

        # With injection at the left dome and sealed caprock, CO2 fills traps in spillpoint order
        # (at minimum the injected trap should fill and spill)
        fill_times = [e.timestamp for e in seq]
        @test all(isfinite, fill_times)
    end

end
