using DGGS
using Test
using GeoDataFrames
using YAXArrays

@testset "DGGS.jl" begin
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "dggs_type" => "ISEA4H",
        "dggs_res_spec" => 2,
        "clip_subset_type" => "WHOLE_EARTH",
        "point_output_type" => "TEXT",
        "point_output_file_name" => "centers"
    )

    d = call_dggrid(meta)
    @test isfile("$(d)/centers.txt")

    @test_throws DomainError Grid("Foo")
    @test_throws DomainError Grid("ISEA", 100, "HEXAGON", 5)

    grid = create_toy_grid()
    @test length(grid.data.data) == 642
    export_cell_boundaries(grid)
    @test GeoDataFrames.read("boundaries.geojson") |> size == (642, 2)
    export_cell_centers(grid)
    @test GeoDataFrames.read("centers.geojson") |> size == (642, 2)
    @test get_cell_ids(grid, 58, 11) == 1
    @test get_cell_ids(grid, 59, 11) == 1
    @test get_geo_coords(grid, 1) == (58.2825256, 11.25)

    geo_cube = Cube("tos_O1_2001-2002.nc")
    cell_cube = get_cell_cube(grid, geo_cube, "Y", "X")
    @test cell_cube.cell_id |> length == 642
    geo_cube2 = get_geo_cube(grid, cell_cube)
    @test isdefined(geo_cube2, :data)

    grid2 = Grid("ISEA4H")
    @test grid2.spec.projection == "ISEA"
    @test grid2.spec.type == "ISEA4H"
    @test grid2.spec.aperture == 4
    @test grid2.spec.resolution == 9
    @test grid2.spec.projection == "ISEA"
    @test grid2.spec.topology == "HEXAGON"

    grid3 = Grid("ISEA", 4, "HEXAGON", 3)
    @test length(grid3.data.data) == 642

    @test get_cell_ids(grid3, 0, 0) == 157
    @test get_cell_ids(grid3, 80, 170) == 313
    @test_throws DomainError get_cell_ids(grid3, 180, 0)
    @test_throws DomainError get_cell_ids(grid3, 0, 200)
    @test get_cell_ids(grid3, -90:5:90, -180:5:180) |> length == 2701

    @test grid3 |> get_cell_boundaries |> size == (642, 2)
    @test grid3 |> get_cell_centers |> size == (642, 2)
end