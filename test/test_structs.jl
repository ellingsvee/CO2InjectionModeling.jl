@testset "Data structures" begin
    # Check that dx, dy, dz are correctly calculated
    domain = Domain3D(10, 10, 5, 100.0, 100.0, 50.0, 800.0, 850.0)
    @test domain.dx ≈ 10.0
    @test domain.dy ≈ 10.0
    @test domain.dz ≈ 10.0

    # Check that the leakage height is computed correctly when shale_pressure_threshold is finite
    rp = ReservoirProperties(0.3, 0.2, 0.1, 50_000.0, 5.0)
    @test rp.sand_porosity == 0.3
    @test isfinite(rp.leakage_height)
    @test rp.leakage_height > 0

    # ...and when shale_pressure_threshold is inf.
    rp_sealed = ReservoirProperties(0.3, 0.2, 0.1, Inf, 5.0)
    @test rp_sealed.leakage_height == Inf

    # Leakage height changes with different CO2 density
    rp_light = ReservoirProperties(0.3, 0.2, 0.1, 50_000.0, 5.0; co2_density=300.0)
    rp_heavy = ReservoirProperties(0.3, 0.2, 0.1, 50_000.0, 5.0; co2_density=600.0)
    # Lighter CO2 → larger density difference → shorter column needed → smaller leakage_height
    @test rp_light.leakage_height < rp_heavy.leakage_height
    # Both should be positive and finite
    @test rp_light.leakage_height > 0
    @test isfinite(rp_heavy.leakage_height)
    # Check formula: h = P / ((ρ_brine - ρ_co2) * g)
    g = 9.81
    @test rp_light.leakage_height ≈ 50_000.0 / ((1020.0 - 300.0) * g)
    @test rp_heavy.leakage_height ≈ 50_000.0 / ((1020.0 - 600.0) * g)
end
