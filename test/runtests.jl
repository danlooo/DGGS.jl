using DGGS
using Test
using GeoDataFrames

@testset "DGGS.jl" begin
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "dggs_type" => "ISEA4H",
        "dggs_res_spec" => 2,
        "clip_subset_type" => "WHOLE_EARTH",
        "point_output_type" => "TEXT",
        "point_output_file_name" => "centers"
    )

    d = dg_call(meta)
    @test isfile("$(d)/centers.txt")

    @test_throws DomainError Grid("Foo")
    @test_throws DomainError Grid("ISEA", 100, "HEXAGON", 5)

    grid = create_toy_grid()
    @test length(grid.data.data) == 642
    export_cell_boundaries(grid)
    @test GeoDataFrames.read("boundaries.geojson") |> size == (642, 2)
    export_cell_centers(grid)
    @test GeoDataFrames.read("centers.geojson") |> size == (642, 1)

    grid2 = Grid("ISEA4H")
    @test grid2.spec.projection == "ISEA"
    @test grid2.spec.type == "ISEA4H"
    @test grid2.spec.aperture == 4
    @test grid2.spec.resolution == 9
    @test grid2.spec.projection == "ISEA"
    @test grid2.spec.topology == "HEXAGON"

    grid3 = Grid("ISEA", 4, "HEXAGON", 3)
    @test length(grid3.data.data) == 642
    @test get_cell_name(grid3, 0, 0) == 157
    @test get_cell_name(grid3, 80, 170) == 289
    @test_throws DomainError get_cell_name(grid3, 180, 0)
    @test_throws DomainError get_cell_name(grid3, 0, 200)

    @test length(get_geo_coords(grid3, [1, 2, 10])) == 3

    @test grid3.spec |> get_cell_boundaries |> size == (642, 2)
    @test length(get_cell_centers(grid3.spec).data) == 642
end