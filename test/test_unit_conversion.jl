@testset "Volume conversion round-trip" begin
    domain = Domain3D(10, 10, 5, 100.0, 100.0, 50.0, 800.0, 850.0)
    rp = ReservoirProperties(0.3, 0.2, 0.1, 50000.0, 5.0)

    vol_physical = 1000.0
    vol_swim = physical_volume_to_swim_volume(vol_physical, rp, domain)
    vol_back = swim_volume_to_physical_volume(vol_swim, rp, domain)
    @test vol_back ≈ vol_physical
end

@testset "Mass conversion" begin
    domain = Domain3D(10, 10, 5, 100.0, 100.0, 50.0, 800.0, 850.0)
    rp = ReservoirProperties(0.3, 0.2, 0.1, 50000.0, 5.0; co2_density=500.0)

    # physical_volume_to_mass: simple multiplication
    @test physical_volume_to_mass(1.0, 500.0) ≈ 500.0
    @test physical_volume_to_mass(2.0, 460.0) ≈ 920.0
    @test physical_volume_to_mass(0.0, 500.0) ≈ 0.0

    # swim_volume_to_mass should equal swim_volume_to_physical_volume * co2_density
    vol_swim = 100.0
    mass = swim_volume_to_mass(vol_swim, rp, domain)
    physical_vol = swim_volume_to_physical_volume(vol_swim, rp, domain)
    @test mass ≈ physical_vol * rp.co2_density

    # Array input
    vols = [10.0, 20.0, 30.0]
    masses = physical_volume_to_mass(vols, 460.0)
    @test masses ≈ vols .* 460.0
end