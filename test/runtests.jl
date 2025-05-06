using DGGS
using Test
using Distances
using ArchGDAL
using YAXArrays
using DimensionalData

@testset "DGGS.jl" begin

    @testset "Cells" begin
        @test Cell(1, 2, 3, 4) isa Cell
        @test Cell{Int32}(1, 2, 3, 4) isa Cell{Int32}

        @test Cell(0, 0, 0, 3) < Cell(1, 0, 0, 3)
        @test_throws ErrorException Cell(0, 0, 0, 3) < Cell(1, 0, 0, 4)

        @test to_geo(1, 1, 1, 8) == to_geo(Cell(1, 1, 1, 8))
    end

    @testset "DGGSArray" begin
        resolution = 3
        i_dim = Dim{:dggs_i}(0:2*2^resolution-1)
        j_dim = Dim{:dggs_j}(0:2^resolution-1)
        n_dim = Dim{:dggs_n}(0:4)
        time_dim = Ti(1:10)
        dim_array = rand(i_dim, j_dim, n_dim, time_dim)
        yax_array = YAXArray(dim_array.dims, dim_array.data)

        @test DGGSArray(dim_array, resolution, "ISEA4D.Penta") isa DGGSArray
        @test DGGSArray(yax_array, resolution, "ISEA4D.Penta") isa DGGSArray
    end

    @testset "Coordinate transformations" begin
        resolution = 20
        geo_points = [(lon, lat) for lat in -90:5:90 for lon in -180:5:180]
        cell_ids = map(x -> to_cell(x..., resolution), geo_points)
        geo_points2 = to_geo.(cell_ids)

        authalic_haversine = Haversine(6371007.18091875)
        dists = colwise(authalic_haversine, geo_points, geo_points2)
        # 99% of points must be < 10m after re-projection
        @test sum(dists .< 10) / length(dists) >= 0.99

        # cell ids must be in bounds
        @test all(map(x -> x.i in 0:2*2^resolution-1, cell_ids))
        @test all(map(x -> x.j in 0:2^resolution-1, cell_ids))
        @test all(map(x -> x.n in 0:4, cell_ids))
    end

    @testset "Integer index" begin
        resolution = 5
        cells = [Cell(i, j, n, resolution) for n in 0:4 for j in 0:2^resolution-1 for i in 0:2*2^resolution-1]
        cells_int = Int64.(cells)
        cells2 = Cell.(cells_int, resolution)

        @test length(cells) == length(cells_int |> unique)
        @test cells_int == 0:length(cells)-1
        @test cells == cells2
    end

    @testset "Convert geo to DGGS" begin
        resolution = 6
        lon_range = X(180:-1:-180)
        lat_range = Y(90:-1:-90)
        geo_data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
        metadata = Dict("standard_name" => "air_temperature", "units" => "K", "description" => "random test data")
        geo_array = YAXArray((lon_range, lat_range), geo_data, metadata)
        dggs_array = to_dggs_array(geo_array, resolution; lon_name=:X, lat_name=:Y)
        geo_array2 = to_geo_array(dggs_array, geo_array.X, geo_array.Y)
        geo_diffs = abs.(geo_array .- geo_array2)

        @test size(geo_array) == size(geo_array2)
        @test all(geo_diffs .< 1.7) # max global deviation
        @test sum(geo_diffs .< 0.1) / length(geo_diffs) >= 0.95
        # alternative methods
        lon_range = -180:180
        lat_range = -90:90
        geo_array3 = to_geo_array(dggs_array, lon_range, lat_range)
        @test size(geo_array3) == (length(lon_range), length(lat_range))
    end
end
