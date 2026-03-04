import SurfaceWaterIntegratedModeling as SWIM

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

const LAYER_MASS_ATOL = 1e-8

"""
Integrate a Vector{WeatherEvent} from t=0 to `t` and return the total CO2
in SWIM units entering the domain. WeatherEvents must have Matrix rain_rate.
"""
function total_received(weather_events::Vector{WeatherEvent}, t::Float64)::Float64
    return compute_total_injected(weather_events, t)
end


@testset "Multi-layer leakage weather events" begin

    # ── Layer 1 simulation setup (leaky caprock) ──────────────────────────────
    layer_idx = 1

    # layers = dome_layers
    # domain = dome_domain
    layers = grf_layers
    domain = grf_domain

    layer1 = layers[layer_idx]
    we1 = convert_injection_event_to_weather_event(
        injection_events[layer_idx], rp_quick, domain)
    seq1, ls1 = fill_sequence_with_leakage(layer1.trap_structure, rp_quick,
        we1; verbose=false)

    # Direct injection for layer 2 (none in this scenario)
    we2_direct = convert_injection_event_to_weather_event(
        injection_events[2], rp_quick, domain)


    @testset "No leakage -> returns direct events unchanged" begin
        # Sealed caprock: nothing drains, result should be copy of direct events
        we1_sealed = convert_injection_event_to_weather_event(
            injection_events[layer_idx], rp_sealed, domain)
        seq_sealed, ls_sealed = fill_sequence_with_leakage(
            layer1.trap_structure, rp_sealed, we1_sealed; verbose=false)

        result = generate_leakage_weather_events(
            seq_sealed, ls_sealed, layer1.trap_structure,
            rp_sealed, rp_sealed, we2_direct)

        # Should be identical to the direct events (nothing leaks)
        @test length(result) == length(we2_direct)
        for (r, d) in zip(result, we2_direct)
            @test r.timestamp == d.timestamp
            if r.rain_rate isa Matrix && d.rain_rate isa Matrix
                @test r.rain_rate ≈ d.rain_rate atol = LAYER_MASS_ATOL
            else
                @test r.rain_rate == d.rain_rate
            end
        end
    end


    @testset "Leaking layer -> non-zero events generated" begin
        @test length(ls1.leakage_records) > 0  # leakage must have occurred

        we2_combined = generate_leakage_weather_events(
            seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

        @test !isempty(we2_combined)

        # Must include at least one event with a non-zero rain rate
        # (at the leakage location)
        any_nonzero = any(any(we.rain_rate .> 0) for we in we2_combined)
        @test any_nonzero

        # Events must be sorted by timestamp
        ts = [we.timestamp for we in we2_combined]
        @test issorted(ts)

        # First timestamp must be ≤ first timestamp of layer-1 sequence
        @test we2_combined[1].timestamp <= seq1[1].timestamp + LAYER_MASS_ATOL
    end


    @testset "Mass conservation: layer-1 output = layer-2 input" begin
        @test length(ls1.leakage_records) > 0

        we2_combined = generate_leakage_weather_events(
            seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

        # Pick timepoints: sim start, mid-drainage, end-drainage, after injection
        t_sim_start = seq1[1].timestamp
        t_leak_start = ls1.leakage_start_time[findfirst(ls1.leaking)]
        t_drain_end = t_leak_start + rp_quick.residual_leakage_time

        timepoints = sort(unique([
            t_sim_start,
            t_leak_start,
            t_leak_start + 0.5 * rp_quick.residual_leakage_time,
            t_drain_end,
            t_drain_end + rp_quick.residual_leakage_time,   # well past drainage
        ]))

        snaps1 = generate_layer_snapshots(
            layer1, layer_idx, seq1, ls1, we1, timepoints)

        for (snap, t) in zip(snaps1, timepoints)
            # CO2 that left layer 1 by time t
            left_layer1 = total_to_next_layer(snap)

            # CO2 received by layer 2 by time t (from combined weather events)
            received_layer2 = total_received(we2_combined, t)

            @test left_layer1 ≈ received_layer2 atol = LAYER_MASS_ATOL
        end
    end


    @testset "Drainage component matches residual drainage" begin
        we2_combined = generate_leakage_weather_events(
            seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

        t_sim_start = seq1[1].timestamp
        t_leak_start = ls1.leakage_start_time[findfirst(ls1.leaking)]
        t_drain_end = t_leak_start + rp_quick.residual_leakage_time

        # After drainage is complete, the drainage component of we2_combined
        # must equal the total residual drained from layer 1.
        t_after = t_drain_end + 1.0

        total_drained_layer1 = compute_total_drained(t_after, ls1)

        # From the combined events, subtract the direct-injection contribution
        # to isolate the leakage portion.
        direct_total = total_received(we2_direct, t_after)
        leakage_total = total_received(we2_combined, t_after) - direct_total

        # leakage_total = total_drained + total_passthrough (both from layer 1)
        # leakage_total ≥ total_drained (passthrough ≥ 0)
        @test leakage_total >= total_drained_layer1 - LAYER_MASS_ATOL

        # And from the layer-1 snapshot, total_to_next_layer = leakage_total
        snaps1 = generate_layer_snapshots(layer1, layer_idx, seq1, ls1, we1, [t_after])
        @test total_to_next_layer(snaps1[1]) ≈ leakage_total atol = LAYER_MASS_ATOL
    end


    @testset "End-to-end: layer-2 simulation runs with generated events" begin
        we2_combined = generate_leakage_weather_events(
            seq1, ls1, layer1.trap_structure, rp_quick, rp_quick, we2_direct)

        layer2 = dome_layers[2]

        # Layer 2 can run with the generated weather events
        seq2, ls2 = fill_sequence_with_leakage(
            layer2.trap_structure, rp_quick, we2_combined; verbose=false)

        @test length(seq2) > 0

        # All timestamps in layer-2 sequence must be finite
        @test all(isfinite(se.timestamp) for se in seq2)

        # Mass conservation within layer 2 at a late timepoint
        t_check = seq2[end].timestamp
        snap2 = generate_layer_snapshot(
            layer2, 2, seq2, ls2, we2_combined, t_check,
            SWIM.trap_states_at_timepoints(layer2.trap_structure, seq2, [t_check]; verbose=false)[1])

        @test snap2.total_stored >= -LAYER_MASS_ATOL
        @test snap2.total_drained >= -LAYER_MASS_ATOL
        @test total_to_next_layer(snap2) >= -LAYER_MASS_ATOL
        @test snap2.total_stored + total_to_next_layer(snap2) ≈ snap2.total_injected atol = LAYER_MASS_ATOL
    end


    @testset "Multiple injection points: independent plumes" begin
        # Use the GRF topology which has multiple traps.
        # Check that generate_leakage_weather_events handles multiple
        # simultaneously leaking traps correctly.
        layer1_grf = grf_layers[layer_idx]
        we1_grf = convert_injection_event_to_weather_event(
            injection_events[layer_idx], rp_quick, grf_domain)

        seq_grf, ls_grf = fill_sequence_with_leakage(
            layer1_grf.trap_structure, rp_quick, we1_grf; verbose=false)

        we2_grf_direct = convert_injection_event_to_weather_event(
            injection_events[2], rp_quick, grf_domain)

        we2_grf = generate_leakage_weather_events(
            seq_grf, ls_grf, layer1_grf.trap_structure,
            rp_quick, rp_quick, we2_grf_direct)

        @test !isempty(we2_grf)
        @test issorted([we.timestamp for we in we2_grf])

        # Mass conservation at a late timepoint
        t_late = seq_grf[end].timestamp + rp_quick.residual_leakage_time
        snap1_grf = generate_layer_snapshot(
            layer1_grf, layer_idx, seq_grf, ls_grf, we1_grf, t_late,
            SWIM.trap_states_at_timepoints(layer1_grf.trap_structure, seq_grf, [t_late]; verbose=false)[1])

        received = total_received(we2_grf, t_late) - total_received(we2_grf_direct, t_late)
        @test received ≈ total_to_next_layer(snap1_grf) atol = LAYER_MASS_ATOL
    end

end
