using DGGS
using DimensionalData
using YAXArrays
using Test

@testset verbose = true "DGGS.jl" begin
    lon_range = X(-180:180)
    lat_range = Y(-90:90)
    level = 6
    data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]

    raster = DimArray(data, (lon_range, lat_range))
    cell_cube = to_cell_cube(raster, level)
    dggs = GridSystem(data, lon_range, lat_range, level)
    dggs = GridSystem(data, -180:180, -90:90, level)

    geo_arr = YAXArray((lon_range, lat_range), data)
    cell_cube2 = to_cell_cube(geo_arr, 5)
    dggs2 = GridSystem(cell_cube2)


    @test typeof(cell_cube) == CellCube
    @test typeof(dggs) == GridSystem
    @test keys(dggs.data) |> maximum == 6

    # convert back and forth
    raster2 = raster |> x -> to_cell_cube(x, level) |> to_geo_cube
    @test name.(dims(raster2)) == (:lon, :lat)

    dggs_path = tempname()
    write(dggs_path, dggs)
    dggs2 = GridSystem(dggs_path)
    rm(dggs_path, recursive=true)


    time_range = 1:12
    level = 6
    bands = 1:3
    data = [band * exp(cosd(lon)) + time * (lat / 90)
            for lat in lat_range, lon in lon_range, time in time_range, band in bands]
    axlist = (
        lat_range,
        lon_range,
        Dim{:time}(time_range),
        Dim{:band}(bands)
    )
    dggs2 = YAXArray(axlist, data) |> x -> to_cell_cube(x, level) |> GridSystem
    @test dggs2 |> CellCube |> x -> x.data |> dims |> length == 5
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
    plot(cell_cube; resolution=100, type=:native)
    plot(dggs; resolution=100)
    plot(dggs2; resolution=100)
    plot(cell_cube, BBox(10, 20, 10, 20); resolution=100)
    plot(dggs, BBox(10, 20, 10, 20); resolution=100)
    plot(dggs2, BBox(10, 20, 10, 20); resolution=100)

    dggs3 = GridSystem("https://s3.bgc-jena.mpg.de:9000/dggs/modis")
    @test typeof(dggs3) == GridSystem
end