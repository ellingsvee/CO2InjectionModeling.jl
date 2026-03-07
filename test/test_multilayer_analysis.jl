function run_multilayer_analysis_tests(scenario::TestScenario)
    @testset "MultiLayerSnapshot: mass conservation and consistency" begin

        layers = scenario.layers
        n_layers = length(layers)

        seqs, leakage_states, weather_events_per_layer =
            fill_layers(layers, scenario.domain, rp_quick, scenario.injection_events)

        # ── fill_layers shape ────────────────────────────────────────────────
        @testset "fill_layers shape" begin
            @test length(seqs) == n_layers
            @test length(leakage_states) == n_layers
            @test length(weather_events_per_layer) == n_layers
        end

        # Build timepoints spanning before-leak → during-drain → after-drain
        # Use the earliest leakage start time across all layers that have leakage.
        leaking_layers = filter(i -> any(leakage_states[i].leaking), 1:n_layers)
        @test !isempty(leaking_layers)  # rp_quick should trigger leakage

        t_sim_start = seqs[1][1].timestamp
        t_leak_start = minimum(
            leakage_states[i].leakage_start_time[findfirst(leakage_states[i].leaking)]
            for i in leaking_layers
        )
        t_drain_end = t_leak_start + rp_quick.residual_leakage_time

        timepoints = sort(unique(filter(isfinite, [
            t_sim_start,
            t_leak_start - 0.1 * rp_quick.residual_leakage_time,
            t_leak_start,
            t_leak_start + 0.25 * rp_quick.residual_leakage_time,
            t_leak_start + 0.5 * rp_quick.residual_leakage_time,
            t_drain_end,
            t_drain_end + 0.5 * rp_quick.residual_leakage_time,
            t_drain_end + rp_quick.residual_leakage_time,
        ])))
        timepoints = filter(t -> t >= t_sim_start, timepoints)

        multi_snaps = generate_multi_layer_snapshots(
            layers, seqs, leakage_states, weather_events_per_layer, timepoints)

        @test length(multi_snaps) == length(timepoints)

        # ── Per-layer mass conservation ──────────────────────────────────────
        @testset "Per-layer conservation" begin
            for msnap in multi_snaps
                for layer_snap in msnap.layers
                    assert_mass_conservation(layer_snap)
                end
            end
        end

        # ── Cross-layer coupling: leakage from layer k = injection into layer k+1 ──
        @testset "Cross-layer coupling" begin
            for msnap in multi_snaps
                for k in 1:(n_layers-1)
                    upward = total_upward_leakage(msnap.layers[k])
                    received = msnap.layers[k+1].total_injected
                    @test upward ≈ received atol = MASS_ATOL
                end
            end
        end

        # ── Snapshot field consistency ────────────────────────────────────────
        @testset "Snapshot field consistency" begin
            for msnap in multi_snaps
                @test msnap.total_injected ≈ msnap.layers[1].total_injected atol = MASS_ATOL
                @test msnap.total_stored ≈ sum(s.total_stored for s in msnap.layers) atol = MASS_ATOL
                @test msnap.total_surface_leakage ≈ total_to_next_layer(msnap.layers[end]) atol = MASS_ATOL
                @test msnap.timestamp == msnap.layers[1].timestamp
                for (i, s) in enumerate(msnap.layers)
                    @test s.layer_idx == i
                    @test s.timestamp == msnap.timestamp
                end
            end
        end

        # ── Monotonicity ──────────────────────────────────────────────────────
        @testset "Monotonicity" begin
            for i in 2:length(multi_snaps)
                @test multi_snaps[i].total_injected >= multi_snaps[i-1].total_injected - MASS_ATOL
                @test multi_snaps[i].total_surface_leakage >= multi_snaps[i-1].total_surface_leakage - MASS_ATOL
            end
        end

        # ── print_summary ─────────────────────────────────────────────────────
        @testset "print_summary" begin
            buf = IOBuffer()
            print_summary(buf, multi_snaps[end])
            output = String(take!(buf))
            @test !isempty(output)
            @test occursin("TOTAL", output)
        end

    end
end
