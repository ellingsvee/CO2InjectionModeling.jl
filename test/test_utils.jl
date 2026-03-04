@testset "Volume to height conversion" begin
    tstruct = DOME_SCENARIO.layers[2].trap_structure

    # Check that the min_topography_elevation is retrieved correctly (single dome: min is at dome center)
    min_topography_elevation = get_min_topography_elevation(1, tstruct)
    @test isfinite(min_topography_elevation)

    # For the single-dome surface, trap 1 has no children, so bottom elevation equals min topography
    trap_bottom_elevation = get_trap_bottom_elevation(1, tstruct)
    @test trap_bottom_elevation == min_topography_elevation
end

@testset "Parent and descendant relationships (GRF topology)" begin
    tstruct = GRF_SCENARIO.layers[1].trap_structure
    n = numtraps(tstruct)
    @test n >= 2

    # Identify root traps (no parent) and child traps (have a parent)
    root_ids  = [id for id in 1:n if isnothing(parentof(tstruct, id))]
    child_ids = [id for id in 1:n if !isnothing(parentof(tstruct, id))]

    # In a closed-BC GRF setup, there should be at least one child trap
    @test !isempty(child_ids)

    # Every child should list at least one parent via get_all_parents
    for child_id in child_ids
        parents = get_all_parents(tstruct, child_id)
        @test !isempty(parents)
    end

    # Every parent trap's descendants should include its direct children
    parent_ids = [id for id in 1:n if !isempty(subtrapsof(tstruct, id))]
    for parent_id in parent_ids
        direct_children = collect(subtrapsof(tstruct, parent_id))
        descendants     = get_all_descendants(tstruct, parent_id)
        @test all(c in descendants for c in direct_children)

        # Parent's trap bottom elevation >= its min topography elevation
        min_elev    = get_min_topography_elevation(parent_id, tstruct)
        bottom_elev = get_trap_bottom_elevation(parent_id, tstruct)
        @test bottom_elev >= min_elev
    end
end
