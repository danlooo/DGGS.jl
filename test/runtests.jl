using DGGS
using Test

@testset verbose = true "DGGS.jl" begin
    lon_range = -180:180
    lat_range = -90:90
    level = 6
    data = [exp(cosd(lon)) + 3(lat / 90) for lat in lat_range, lon in lon_range]
    dggs = GridSystem(data, lon_range, lat_range, level)

    @test keys(dggs.data) |> maximum == 6
end