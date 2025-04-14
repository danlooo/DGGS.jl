using PentaCube
using Test
using Distances

@testset "PentaCube.jl" begin

    @testset "Cells" begin
        @test Cell(1, 2, 3, 4) isa Cell
        @test Cell{Int32}(1, 2, 3, 4) isa Cell{Int32}

        @test Cell(0, 0, 0, 3) < Cell(1, 0, 0, 3)
        @test_throws ErrorException Cell(0, 0, 0, 3) < Cell(1, 0, 0, 4)

        @test to_geo(1, 1, 1, 8) == to_geo(Cell(1, 1, 1, 8))
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
end
