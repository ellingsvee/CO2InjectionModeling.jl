# @testset "Get parents and descendants" begin
#     pass
# end

@testset "Volume to height conversion" begin
    tstruct = layers[1].trap_structure

    # Check that the min_topography_elevation is retrieved correctly
    min_topography_elevation = get_min_topography_elevation(1, tstruct)
    @test min_topography_elevation == minimum(_s1)

    # For this surface, the trap bottom elevation should be the same as the min topography elevation since there are no child traps
    trap_bottom_elevation = get_trap_bottom_elevation(1, tstruct)
    @test trap_bottom_elevation == min_topography_elevation

end