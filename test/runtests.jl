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
    end

    @testset "Import NetCDF" begin
        using YAXArrays
        using NetCDF
        using Downloads

        url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
        file_path = tempname()
        filename = Downloads.download(url, file_path)
        geo_cube_raw = YAXArrays.Cube(file_path)

        data = circshift(geo_cube_raw[:, :, 1], 90)
        latitudes = geo_cube_raw.Y
        longitudes = geo_cube_raw.X .- 180
        geo_cube = GeoCube(data, latitudes, longitudes)
        @test size(geo_cube.data) == (180, 170)
        dggs = DgGlobalGridSystem(geo_cube, 3, :isea, 4, :hexagon)
        @test length(dggs[1]) == 12
        rm(file_path)
    end
end