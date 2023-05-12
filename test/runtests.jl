using DGGS
using Test

@testset "DGGS.jl" begin
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "dggs_type" => "ISEA4H",
        "clip_subset_type" => "WHOLE_EARTH",
        "cell_output_type" => "GEOJSON",
        "cell_output_file_name" => "out.geojson"
    )

    d = dg_call(meta)

    grid = Grid(DGGS.ISEA4H)
    @test grid.projection == DGGS.ISEA
    @test grid.type == "ISEA4H"
    @test grid.aperture == 4
    @test grid.resolution == 9
    @test grid.projection == DGGS.ISEA
    @test grid.topology == DGGS.HEXAGON
    cells = generate_cells(grid)

    grid2 = Grid(ISEA, 4, HEXAGON, 9)
    @test grid == grid22
end