using DGGS
using Test
using GeoJSON

@testset "DGGS.jl" begin
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "dggs_type" => "ISEA4H",
        "dggs_res_spec" => 2,
        "clip_subset_type" => "WHOLE_EARTH",
        "cell_output_type" => "GEOJSON",
        "cell_output_file_name" => "out"
    )

    d = dg_call(meta)
    @test length(GeoJSON.read(open("$(d)/out.geojson"), ndim=3)) == 162

    grid = Grid("ISEA4H")
    @test grid.projection == "ISEA"
    @test grid.type == "ISEA4H"
    @test grid.aperture == 4
    @test grid.resolution == 9
    @test grid.projection == "ISEA"
    @test grid.topology == "HEXAGON"

    @test_throws DomainError Grid("Foo")
    @test_throws DomainError Grid("ISEA", 100, "HEXAGON", 5)
end