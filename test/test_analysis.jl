using SurfaceWaterIntegratedModeling: numtraps

# Tolerance used for floating-point mass balance comparisons.
const MASS_BALANCE_ATOL = 1e-8


"""
Expected total CO2 injected into a layer up to time `t` (in SWIM units),
computed directly from the injection schedule constants and unit conversion.

This is independent of the WeatherEvent / SpillEvent machinery and serves as
a ground-truth reference for `snap.total_injected`.
"""
function expected_total_injected_swim(
    t::Float64,
    rp::ReservoirProperties,
    domain::Domain3D
)::Float64
    active_duration = clamp(t - INJECTION_START_T, 0.0, INJECTION_END_T - INJECTION_START_T)
    physical_injected = TOTAL_INJECTION_RATE * active_duration
    return physical_volume_to_swim_volume(physical_injected, rp, domain)
end

"""
Assert that a LayerSnapshot satisfies exact mass conservation.
"""
function assert_mass_conservation(snap::LayerSnapshot; atol=MASS_BALANCE_ATOL)
    @test snap.total_stored >= -atol
    @test snap.total_drained >= -atol
    @test total_to_next_layer(snap) >= -atol
    @test total_passthrough(snap) >= -atol
    # True by definition, but guards against implementation errors
    # where total_stored or total_injected are computed inconsistently.
    @test snap.total_stored + total_to_next_layer(snap) ≈ snap.total_injected atol = atol
end

@testset "Layer snapshots and mass conservation" begin

    # layers = dome_layers
    # domain = dome_domain
    layers = grf_layers
    domain = grf_domain

    @testset "Sealed caprock – no leakage" begin
        layer_idx = 1
        layer = layers[layer_idx]
        we_events = convert_injection_event_to_weather_event(
            injection_events[layer_idx], rp_sealed, domain)

        seq, leakage_state = fill_sequence_with_leakage(layer.trap_structure, rp_sealed,
            we_events; verbose=false)

        # Injection runs from t=0 to t=10, then rate drops to 0.
        # Snapshot at t=5 (mid-injection) and t=20 (well after injection ends).
        t_sim_start = seq[1].timestamp
        timepoints = [t_sim_start + 5.0, t_sim_start + 20.0]

        snaps = generate_layer_snapshots(layer, layer_idx, seq, leakage_state,
            we_events, timepoints)

        for snap in snaps
            assert_mass_conservation(snap)
        end

        # With sealed caprock nothing should ever drain or leave the layer
        for snap in snaps
            @test snap.total_drained ≈ 0.0 atol = MASS_BALANCE_ATOL
            @test total_to_next_layer(snap) ≈ 0.0 atol = MASS_BALANCE_ATOL
            @test snap.total_stored ≈ snap.total_injected atol = MASS_BALANCE_ATOL
        end

        # CO2 should accumulate monotonically while injection is active
        @test snaps[1].total_injected <= snaps[2].total_injected + MASS_BALANCE_ATOL

        # Verify total_injected matches the known physical injection schedule
        for snap in snaps
            expected = expected_total_injected_swim(snap.timestamp, rp_sealed, domain)
            @test snap.total_injected ≈ expected atol = MASS_BALANCE_ATOL
        end
    end


    @testset "Leaking layer: Mass conservation with drainage" begin
        layer_idx = 1
        layer = layers[layer_idx]
        we_events = convert_injection_event_to_weather_event(
            injection_events[layer_idx], rp_quick, domain)

        seq, leakage_state = fill_sequence_with_leakage(layer.trap_structure, rp_quick,
            we_events; verbose=false)

        @test length(leakage_state.leakage_records) > 0  # leakage must have occurred

        t_sim_start = seq[1].timestamp
        # Sample before, during, and after the drainage window
        t_leak_start = leakage_state.leakage_start_time[findfirst(leakage_state.leaking)]
        timepoints = sort(unique([
            t_sim_start,
            t_leak_start + 0.5 * rp_quick.residual_leakage_time,
            t_leak_start + rp_quick.residual_leakage_time,           # drainage complete
            t_leak_start + 2.0 * rp_quick.residual_leakage_time,
        ]))

        snaps = generate_layer_snapshots(layer, layer_idx, seq, leakage_state,
            we_events, timepoints)

        for snap in snaps
            assert_mass_conservation(snap)
        end

        # total_to_next_layer should be non-decreasing over time (CO2 doesn't flow back)
        for i in 2:length(snaps)
            @test total_to_next_layer(snaps[i]) >= total_to_next_layer(snaps[i-1]) - MASS_BALANCE_ATOL
        end

        # total_drained should also be non-decreasing (drainage only goes forward)
        for i in 2:length(snaps)
            @test snaps[i].total_drained >= snaps[i-1].total_drained - MASS_BALANCE_ATOL
        end

        # Verify total_injected matches the known physical injection schedule
        for snap in snaps
            expected = expected_total_injected_swim(snap.timestamp, rp_quick, domain)
            @test snap.total_injected ≈ expected atol = MASS_BALANCE_ATOL
        end

        # After the drainage window the residual drainage should be complete
        snap_after = snaps[end]
        for trap in 1:numtraps(layer.trap_structure)
            if leakage_state.draining[trap]
                t_drain_end = leakage_state.leakage_start_time[trap] + rp_quick.residual_leakage_time
                if snap_after.timestamp >= t_drain_end
                    expected_final_vol = leakage_state.initial_volume_at_leak[trap] *
                                         rp_quick.sand_residual_co2_saturation
                    @test snap_after.trap_volumes[trap] ≈ expected_final_vol atol = MASS_BALANCE_ATOL
                end
            end
        end
    end


    @testset "Snapshot fields are self-consistent" begin
        layer_idx = 1
        layer = layers[layer_idx]
        we_events = convert_injection_event_to_weather_event(
            injection_events[layer_idx], rp_quick, domain)
        seq, leakage_state = fill_sequence_with_leakage(layer.trap_structure, rp_quick,
            we_events; verbose=false)

        t = seq[1].timestamp + 8.0
        snap = generate_layer_snapshot(layer, layer_idx, seq, leakage_state, we_events, t, trap_states_at_timepoints(layer.trap_structure, seq, [t]; verbose=false)[1])

        # total_stored must equal sum of per-trap volumes
        @test snap.total_stored ≈ sum(snap.trap_volumes) atol = MASS_BALANCE_ATOL

        # trap_volumes must be non-negative
        @test all(snap.trap_volumes .>= -MASS_BALANCE_ATOL)

        # leaking traps must also be draining (leaking implies draining)
        for trap in eachindex(snap.trap_leaking)
            if snap.trap_leaking[trap]
                @test snap.trap_draining[trap]
            end
        end

        # Snapshot name and index are correctly set
        @test snap.layer_idx == layer_idx
        @test snap.layer_name == layer.name
        @test snap.timestamp == t
    end

end
