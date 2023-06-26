using DGGS
using Test

@testset verbose = true "DGGS.jl" begin
    @testset "DGGRID" begin
        meta = Dict(
            "dggrid_operation" => "GENERATE_GRID",
            "dggs_type" => "ISEA4H",
            "dggs_res_spec" => 2,
            "clip_subset_type" => "WHOLE_EARTH",
            "point_output_type" => "TEXT",
            "point_output_file_name" => "centers"
        )

        d = DGGS.call_dggrid(meta)
        @test isfile("$(d)/centers.txt")

        @test DGGS.get_dggrid_grid_table(:isea4h, 3) |> size == (642, 3)
    end

    @testset "Grids" begin
        @test_throws Exception Grid(:isea, 100, :hexagon, 5)

        grid = create_toy_grid()
        @test length(grid) == 642

        @test get_cell_boundaries(grid) |> size == (642, 2)
        @test get_cell_centers(grid) |> size == (642, 2)

        @test get_cell_ids(grid, 58, 11) == 1
        @test get_cell_ids(grid, 59, 11) == 1
        @test get_geo_coords(grid, 1) == (58.2825256, 11.25)

        grid2 = Grid([-170 -80; -165.12 81.12; -160 90]')
        @test length(grid2) == 3

        grid3 = DgGrid(:isea, 4, :hexagon, 3)
        @test length(grid3) == 642

        @test get_cell_ids(grid3, 0, 0) == 157
        @test get_cell_ids(grid3, 80, 170) == 313
        @test get_cell_ids(grid3, -90:5:90, -180:5:180) |> length == 2701

        @test grid3 |> get_cell_boundaries |> size == (642, 2)
        @test grid3 |> get_cell_centers |> size == (642, 2)

        grid4 = DgGrid(:superfund, 3)
        @test grid4.projection == :fuller
        @test isnothing(grid4.aperture)
        @test length(grid4) == 12962
    end

    @testset "Cubes" begin
        lon_range = -180:180
        lat_range = -90:90
        data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
        geo_cube = GeoCube(data, lat_range, lon_range)
        @test isdefined(geo_cube, :data)
        @test eltype(geo_cube) <: Real

        grid = create_toy_grid()
        cell_cube = CellCube(geo_cube, grid)
        @test length(cell_cube) == 642
        geo_cube2 = GeoCube(cell_cube)
        @test isdefined(geo_cube2, :data)

        # test non spatial index dimension, e.g. time
        using YAXArrays
        using NetCDF
        cube_3 = Cube("tos_O1_2001-2002.nc")
        geo_cube3 = GeoCube(cube_3)
        @test size(geo_cube3.data) == (180, 170, 24)
        cell_cube3 = CellCube(geo_cube3, grid)
        @test size(cell_cube3.data) == (329, 24)

        # test conversion back to geoCube again
        geo_cube3 = GeoCube(cell_cube3)
        @test size(geo_cube3.data) == size(geo_cube3.data)
    end

    @testset "GridSystems" begin
        lon_range = -180:180
        lat_range = -90:90
        data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
        geo_cube = GeoCube(data, lat_range, lon_range)

        dggs = DgGlobalGridSystem(geo_cube, 3)
        @test length(dggs) == 3
        @test [x for x in dggs] |> length == 3 # test iterator
        @test length(dggs[1]) == 12
        @test length(dggs[1][1:10]) == 10

        dggs2 = DgGlobalGridSystem(geo_cube, :superfund, 3)
        @test length(dggs2[1]) == 42
    end
end