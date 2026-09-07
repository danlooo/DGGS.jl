using DGGS
using Test
using Distances
using ArchGDAL
using YAXArrays
using DimensionalData
using Extents
using Makie
using Zarr
using Statistics

resolution = 5
lon_range = X(180:-1:-180)
lat_range = Y(90:-1:-90)
geo_data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
properties = Dict("standard_name" => "air_temperature", "units" => "K", "description" => "random test data")
geo_array = YAXArray((lon_range, lat_range), geo_data, properties)
geo_ds = Dataset(Gray=geo_array)
dggs_array = to_dggs_array(geo_array, resolution, "EPSG:4326")

properties2 = Dict("standard_name" => "precipitation")
geo_array2 = YAXArray((lon_range, lat_range), geo_data, properties2)
dggs_array2 = to_dggs_array(geo_array2, resolution, "EPSG:4326")

dggs_ds = DGGSDataset(dggs_array, dggs_array2)

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
        i_dim = Dim{:dggs_i}(0:(2*2^resolution-1))
        j_dim = Dim{:dggs_j}(0:(2^resolution-1))
        n_dim = Dim{:dggs_n}(0:4)
        time_dim = Ti(1:10)
        dim_array = rand(i_dim, j_dim, n_dim, time_dim)
        yax_array = YAXArray(dim_array.dims, dim_array.data)
        geo_bbox = Extent(X=(-180, 180), Y=(-90, 90))

        @test DGGSArray(dim_array, resolution, "ISEA4D.Penta", geo_bbox) isa DGGSArray
        @test DGGSArray(yax_array, resolution, "ISEA4D.Penta", geo_bbox) isa DGGSArray
        @test DGGSArray(yax_array, resolution, "ISEA4D.Penta", geo_bbox)[Cell(1, 2, 3, resolution)] isa YAXArray
    end

    @testset "TileArray" begin
        a = TileArray(0, (100, 100), (10, 10))
        @test all(a .== 0)
        a[1, 1] = 1
        @test a[1, 1] == 1
        @test all(a[2:100, :] .== 0)
        a[1, 1] = 0
        @test all(a .== 0)

        a = TileArray{Union{Missing,Int}}(0, (100, 100), (10, 10))
        a[1, 1] = 1
        a[100, 100] = missing
        @test a[1, 1] == 1
        @test ismissing(a[100, 100])

        resolution = 25
        a = DGGSArray(resolution)
        @test a.data isa TileArray
        @test size(a) == (2 * 2^resolution, 2^resolution, 5)
        a[2^resolution, 1, :] = 1:5
        @test a.data[2^resolution, 1, :] == 1:5

        @test YAXArray(TileArray{Union{Missing,Float64}}(0.0, (100, 100), (10, 10))) isa YAXArray
        @test YAXArray(TileArray(0, (100, 100), (10, 10))) isa YAXArray
    end

    @testset "TileArray mapreduce" begin
        # Test mapreduce with non-missing default
        a = TileArray(0, (10, 10), (5, 5))
        @test mapreduce(identity, +, a) == 0
        @test mapreduce(x -> x + 1, +, a) == 100
        a[1, 1] = 42
        @test mapreduce(identity, +, a) == 42
        @test maximum(a) == 42
        @test minimum(a) == 0

        # Test mapreduce with missing default
        b = TileArray{Union{Missing,Int}}(missing, (10, 10, 10), (5, 5, 5))
        b[1, 1, 1] = 2
        b[2, 1, 1] = 3
        @test maximum(b) == 3
        @test minimum(skipmissing(b)) == 2
        @test sum(skipmissing(b)) == 5
        @test maximum(skipmissing(b)) == 3
        @test b |> skipmissing |> maximum == 3

        # Test mapreduce with dims
        c = TileArray(1, (4, 6), (4, 6))
        c[1, 1] = 10
        c[2, 3] = 5
        dense = [c[i, j] for i in 1:4, j in 1:6]
        @test mapreduce(identity, +, c; dims=1) == sum(dense; dims=1)
        @test mapreduce(identity, +, c; dims=2) == sum(dense; dims=2)
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
        @test all(map(x -> x.i in 0:(2*2^resolution-1), cell_ids))
        @test all(map(x -> x.j in 0:(2^resolution-1), cell_ids))
        @test all(map(x -> x.n in 0:4, cell_ids))
    end

    @testset "Integer index" begin
        resolution = 5
        cells = [Cell(i, j, n, resolution) for n in 0:4 for j in 0:(2^resolution-1) for i in 0:(2*2^resolution-1)]
        cells_int = Int64.(cells)
        cells2 = Cell.(cells_int, resolution)

        @test length(cells) == length(cells_int |> unique)
        @test cells_int == 0:(length(cells)-1)
        @test cells == cells2
    end

    @testset "Convert geo to DGGS" begin
        geo_array2 = to_geo_array(dggs_array, geo_array.X, geo_array.Y)
        geo_diffs = abs.(geo_array .- geo_array2)

        @test size(geo_array) == size(geo_array2)
        @test all(geo_diffs .< 2.0) # max global deviation
        @test sum(geo_diffs .< 0.2) / length(geo_diffs) >= 0.95

        # alternative methods
        lon_range = -180:180
        lat_range = -90:90
        geo_array3 = to_geo_array(dggs_array, lon_range, lat_range)
        @test size(geo_array3) == (length(lon_range), length(lat_range))

        # other crs
        geo_array3 = open_dataset("data/geomatrix.tif").Gray
        projection = geo_array3.properties["projection"]
        dggs_array3 = to_dggs_array(geo_array3, 10, projection)
        @test dggs_array3 isa DGGSArray

        # other agg_func
        dggs_array4 = to_dggs_array(geo_array3, 10, projection; agg_func=median)
        @test dggs_array4 isa DGGSArray
    end

    @testset "Plot" begin
        fig = plot(dggs_array)
        cb = filter(x -> x isa Colorbar, fig.content)[1]

        @test fig isa Figure
        @test cb.label[] == "air_temperature"
    end

    @testset "Open and save DGGSArray" begin
        @test dggs_array == dggs_array |> YAXArray |> DGGSArray
        temp_dir = tempname() * ".dggs.zarr"
        save_dggs_array(temp_dir, dggs_array)
        dggs_array2 = open_dggs_array(temp_dir)
        @test dggs_array == dggs_array2
        @test name(dggs_array) == name(dggs_array2)
        rm(temp_dir, recursive=true)
    end


    @testset "Open and save DGGSDataset" begin
        @test dggs_ds == dggs_ds |> Dataset |> DGGSDataset
        temp_dir = tempname() * ".dggs.zarr"
        save_dggs_dataset(temp_dir, dggs_ds)
        dggs_ds2 = open_dggs_dataset(temp_dir)
        @test dggs_ds == dggs_ds2
        rm(temp_dir, recursive=true)
    end


    @testset "DGGSDataset" begin
        resolution = 3
        i_dim = Dim{:dggs_i}(0:(2*2^resolution-1))
        j_dim = Dim{:dggs_j}(0:(2^resolution-1))
        n_dim = Dim{:dggs_n}(0:4)
        time_dim = Ti(1:10)
        dim_array = rand(i_dim, j_dim, n_dim, time_dim)

        a1 = DGGSArray(dim_array, resolution; name=:red)
        a2 = DGGSArray(dim_array, resolution; name=:blue)
        a3 = DGGSArray(dim_array, resolution; name=:green)
        a4 = DGGSArray(rand(i_dim, j_dim, n_dim), resolution; name=:height)
        a5 = DGGSArray(dim_array, resolution)

        ds1 = DGGSDataset(a1)
        ds2 = DGGSDataset(a1, a2, a3, a4)
        ds3 = DGGSDataset(a5)

        @test ds1 isa DGGSDataset
        @test ds2 isa DGGSDataset
        @test ds3 isa DGGSDataset
        @test ds1.red isa DGGSArray
        @test length(keys(ds1)) == 1
        @test length(keys(ds2)) == 4
        @test ds2.resolution == ds2.blue.resolution
        @test ds2.dggsrs == ds2.blue.dggsrs
        @test_throws ErrorException DGGSDataset(a1, a1)
    end

    @testset "DGGSPyramid" begin
        A = [1 1 2 2; 1 1 2 2; 3 3 4 4; 3 3 4 4]
        @test DGGS.coarsen(A, (2, 2)) == [1 2; 3 4]

        dggs_p = to_dggs_pyramid(dggs_ds)
        @test dggs_p isa DGGSPyramid
        @test dggs_p.dggs_s3 isa DGGSDataset
        @test dggs_p.dggs_s3.air_temperature isa DGGSArray
        @test length(dggs_p.branches) == dggs_ds.resolution
        @test dggs_p.dggsrs == dggs_ds.dggsrs
        @test dggs_p.bbox == dggs_ds.bbox
        @test dggs_p.dggs_s3.resolution == dggs_p[3].resolution

        @testset "all values are present" begin
            data = collect(dggs_p[3].precipitation)
            for n in 1:5
                @test length(data[:, :, n]) == length(filter(!ismissing, data[:, :, n]))
            end
        end

        @testset "save and open pyramid" begin
            temp_dir = tempname() * ".dggs.zarr"
            @info temp_dir
            save_dggs_pyramid(temp_dir, dggs_p)
            dggs_p2 = open_dggs_pyramid(temp_dir)
            @test dggs_p.bbox == dggs_p2.bbox
            @test dggs_p.dggsrs == dggs_p2.dggsrs
            @test length(dggs_p.data) == length(dggs_p2.data)
            @test all(keys(dggs_p.data) .== keys(dggs_p2.data))

            # both layers must be present after save and open
            @test name(dggs_p.dggs_s3.air_temperature) == name(dggs_p2.dggs_s3.air_temperature)
            @test name(dggs_p.dggs_s3.precipitation) == name(dggs_p2.dggs_s3.precipitation)
            rm(temp_dir, recursive=true)
        end

        # pyramid from just one array
        dggs_p2 = to_dggs_pyramid(dggs_array)
        @test all([dggs_p2[x] isa DGGSDataset for x in 1:3])
        @test all([dggs_p2[x].air_temperature isa DGGSArray for x in 1:3])

        # pyramid from a subset
        geo_array4 = geo_array[X=20:30, Y=17:22]
        p = to_dggs_pyramid(geo_array4, resolution, "EPSG:4326")

        # plot colormap pyramids
        @test plot(dggs_p, :air_temperature) isa Figure

        # plot RGB pyramids
        dggs_ds_rgb = DGGSDataset(
            rebuild(dggs_array; name=:Red, metadata=Dict()),
            rebuild(dggs_array; name=:Green, metadata=Dict()),
            rebuild(dggs_array; name=:Blue, metadata=Dict()),
        )
        dggs_p_rgb = to_dggs_pyramid(dggs_ds_rgb)
        @test plot(dggs_p_rgb, :Red, :Green, :Blue) isa Figure
    end

    @testset "Init global dataset" begin
        temp_dir = tempname() * ".dggs.zarr"
        global_dggs_ds = DGGS.init_global_dggs_dataset(geo_ds, resolution, temp_dir)
        @test size(global_dggs_ds) == (2 * 2^resolution, 2^resolution, 5)
        rm(temp_dir, recursive=true)
    end
end
