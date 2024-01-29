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

    scene1 = plot(cell_cube)
    scene2 = plot(dggs)

    using YAXArrays
    lon_range = -180:180
    lat_range = -90:90
    time_range = 1:12
    level = 6
    bands = 1:3
    data = [band * exp(cosd(lon)) + time * (lat / 90) for lat in lat_range, lon in lon_range, time in time_range, band in bands]
    axlist = (
        Dim{:lat}(lat_range),
        Dim{:lon}(lon_range),
        Dim{:time}(time_range),
        Dim{:band}(bands)
    )
    dggs2 = YAXArray(axlist, data) |> GeoCube |> x -> CellCube(x, level) |> GridSystem
    @test dggs2[6][band=1, time=2] |> typeof == CellCube
    @test dggs2[6][band=2, time=2].data.axes |> length == 3
end