using DGGS
using Test

@testset verbose = true "DGGS.jl" begin
    lon_range = -180:180
    lat_range = -90:90
    level = 6
    data = [exp(cosd(lon)) + 3(lat / 90) for lat in lat_range, lon in lon_range]

    geo_cube = GeoCube(data, lon_range, lat_range)
    cell_cube = CellCube(geo_cube, level)
    dggs = GridSystem(data, lon_range, lat_range, level)

    @test typeof(geo_cube) == GeoCube
    @test typeof(cell_cube) == CellCube
    @test typeof(dggs) == GridSystem
    @test geo_cube |> x -> CellCube(x, 8) |> GeoCube |> typeof == GeoCube
    @test keys(dggs.data) |> maximum == 6

    dggs_path = tempname()
    write(dggs_path, dggs)
    dggs2 = GridSystem(dggs_path)
    rm(dggs_path, recursive=true)

    using YAXArrays
    lon_range = -180:180
    lat_range = -90:90
    time_range = 1:12
    level = 6
    bands = 1:3
    data = [band * exp(cosd(lon)) + time * (lat / 90)
            for lat in lat_range, lon in lon_range, time in time_range, band in bands]
    axlist = (
        Dim{:lat}(lat_range),
        Dim{:lon}(lon_range),
        Dim{:time}(time_range),
        Dim{:band}(bands)
    )
    dggs2 = YAXArray(axlist, data) |> GeoCube |> x -> CellCube(x, level) |> GridSystem
    @test dggs2[6][band=1, time=2] |> typeof == CellCube
    @test dggs2[6][band=2, time=2].data.axes |> length == 3

    @test transform_points(-180:180, -90:90, 7) |> size == (361, 181)
    @test transform_points(-180, -90, 7) |> size == (1, 1)
    @test transform_points([(-180, -90), (0, 0), (180, 90)], 7) |> length == 3
    @test transform_points([Q2DI(1, 1, 1), Q2DI(2, 2, 2)], 7) |> length == 2

    # converting back and forth must be approx the same at higher levels
    geo_coords = [(lon, lat) for lon in -175:5:175, lat in -85:5:85] |> vec
    level = 12
    fwd_rev_geo_goords = geo_coords |> x -> transform_points(x, level) |> x -> transform_points(x, level)
    @test geo_coords == map(x -> (round(x[1]), round(x[2])), fwd_rev_geo_goords)

    using GLMakie
    plot(cell_cube; resolution=100)
    plot(cell_cube; resolution=100)
    plot(dggs; resolution=100)
    plot(dggs2; resolution=100)
    plot(cell_cube, BBox(10, 20, 10, 20); resolution=100)
    plot(dggs, BBox(10, 20, 10, 20); resolution=100)
    plot(dggs2, BBox(10, 20, 10, 20); resolution=100)
end