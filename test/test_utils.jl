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
    root_ids = [id for id in 1:n if isnothing(parentof(tstruct, id))]
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
        descendants = get_all_descendants(tstruct, parent_id)
        @test all(c in descendants for c in direct_children)

        # Parent's trap bottom elevation >= its min topography elevation
        min_elev = get_min_topography_elevation(parent_id, tstruct)
        bottom_elev = get_trap_bottom_elevation(parent_id, tstruct)
        @test bottom_elev >= min_elev
    end
end

@testset "spread_rate!" begin
    # radius=0: single cell
    m = zeros(5, 5)
    spread_rate!(m, (3, 3), 10.0, 0)
    @test m[3, 3] == 10.0
    @test sum(m) == 10.0

    # radius=1: disc pattern (center + 4 neighbors = 5 cells)
    m = zeros(5, 5)
    spread_rate!(m, (3, 3), 10.0, 1)
    @test sum(m) ≈ 10.0
    @test m[3, 3] ≈ 2.0   # 10/5
    @test m[2, 3] ≈ 2.0
    @test m[4, 3] ≈ 2.0
    @test m[3, 2] ≈ 2.0
    @test m[3, 4] ≈ 2.0
    @test m[2, 2] == 0.0   # diagonal is outside radius=1

    # radius=2: larger disc
    m = zeros(7, 7)
    spread_rate!(m, (4, 4), 13.0, 2)
    @test sum(m) ≈ 13.0
    n_cells = count(x -> x > 0, m)
    @test n_cells == 13  # disc of radius 2 has 13 cells

    # Edge: location near boundary
    m = zeros(5, 5)
    spread_rate!(m, (1, 1), 6.0, 2)
    @test sum(m) ≈ 6.0
    n_cells = count(x -> x > 0, m)
    @test n_cells < 13  # clipped by boundary

    # Accumulation: calling twice adds up
    m = zeros(5, 5)
    spread_rate!(m, (3, 3), 5.0, 0)
    spread_rate!(m, (3, 3), 3.0, 0)
    @test m[3, 3] == 8.0
end

@testset "create_injection_rate with radius" begin
    layers = DOME_SCENARIO.layers
    pad = layers[1].pad_width

    # radius=0 (default): single cell
    rate_mat = create_injection_rate(layers, (10, 10), 100.0)
    @test rate_mat[10+pad, 10+pad] == 100.0
    @test sum(rate_mat) == 100.0

    # radius=1: spread across disc
    rate_mat = create_injection_rate(layers, (10, 10), 100.0; radius=1)
    @test sum(rate_mat) ≈ 100.0
    @test rate_mat[10+pad, 10+pad] ≈ 20.0  # 100/5
end
