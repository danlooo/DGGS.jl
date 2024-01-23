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

    write("example.dggs", dggs)

    scence = plot(cell_cube)
    @test typeof(scence) == Scene
end