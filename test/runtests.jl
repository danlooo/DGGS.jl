using DGGS
using Test
using GeoJSON

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

    grid = Grid("ISEA4H")
    @test grid.spec.projection == "ISEA"
    @test grid.spec.type == "ISEA4H"
    @test grid.spec.aperture == 4
    @test grid.spec.resolution == 9
    @test grid.spec.projection == "ISEA"
    @test grid.spec.topology == "HEXAGON"

    @test_throws DomainError Grid("Foo")
    @test_throws DomainError Grid("ISEA", 100, "HEXAGON", 5)

    grid2 = Grid("ISEA", 4, "HEXAGON", 3)
    @test length(grid2.data.data) == 642
    @test cell_name(grid2, 0, 0) == 157
    @test cell_name(grid2, 80, 170) == 289
    @test_throws DomainError cell_name(grid2, 180, 0)
    @test_throws DomainError cell_name(grid2, 0, 200)
    @test length(geo_coords(grid2, [1, 2, 10])) == 3
end