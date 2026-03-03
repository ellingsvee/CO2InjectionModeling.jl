@testset "Volume conversion round-trip" begin
    domain = Domain3D(10, 10, 5, 100.0, 100.0, 50.0, 800.0, 850.0)
    rp = ReservoirProperties(0.3, 0.2, 0.1, 50000.0, 5.0)

    vol_physical = 1000.0
    vol_swim = physical_volume_to_swim_volume(vol_physical, rp, domain)
    vol_back = swim_volume_to_physical_volume(vol_swim, rp, domain)
    @test vol_back ≈ vol_physical
end