@testset "Fill sequence for layer" begin
    # A super simple test to check that the fill sequence algorithm runs without error for a single layer with a sealed shale layer above it. We don't check the actual output here, just that it runs and produces a non-empty sequence.
    layer_idx = 1
    tstruct = layers[layer_idx].trap_structure
    we_events = convert_injection_event_to_weather_event(injection_events[layer_idx], rp_sealed, domain)
    seq = fill_sequence_with_leakage(tstruct, we_events, verbose=true)
    @test length(seq) > 0

end
