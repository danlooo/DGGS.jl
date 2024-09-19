using DGGS
using DimensionalData
using YAXArrays
using GLMakie
using Test
using Rasters
using Random

Random.seed!(1337)

#
# Ppen arrays, layers, and pyramids
#

p = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/example-ccsm3")
l = p[4]
a = l.tas

@test a.id == :tas
@test length(a.attrs) > length(p.attrs)
@test (setdiff(p.attrs, l.attrs) .|> x -> x.first) == ["dggs_levels", "dggs_level"] # same global attrs expect DGGS level

@test a[10, 1, 1] isa YAXArray
@test (a[10, 1, 1] .== a[Q2DI(10, 1, 1)]) |> collect |> all

#
# Convert lat/lon rasters into a DGGS
#

lon_range = X(-180:180)
lat_range = Y(-90:90)
level = 6
data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
raster = DimArray(data, (lon_range, lat_range))
a2 = to_dggs_array(raster, level; lon_name=:X, lat_name=:Y)
raster2 = to_geo_array(a2, lon_range.val, lat_range.val)
@test a2.level == level
# low loss after converting back
@test all(abs.(Matrix(raster) .- Matrix(raster2)) .< 1.8)

# test aggregation methods
data_bool = [exp(cosd(lon)) + 3(lat / 90) > 1 ? true : missing for lon in lon_range, lat in lat_range]
raster_bool = DimArray(data_bool, (lon_range, lat_range))
@test Float64 in to_dggs_pyramid(raster_bool, level; lon_name=:X, lat_name=:Y, agg_type=:convert)[level-1].layer.data |> eltype |> Base.uniontypes
@test Bool in to_dggs_pyramid(raster_bool, level; lon_name=:X, lat_name=:Y, agg_type=:round)[level-1].layer.data |> eltype |> Base.uniontypes

# load netcdf geo raster
# reformat lon axes from [0,360] to [-180,180]
# skip mask
geo_ds = open_dataset("sresa1b_ncar_ccsm3-example.nc")
geo_ds.axes[:lon] = vcat(range(0, 180; length=128), range(-180, 0; length=128)) |> Dim{:lon}
arrs = Dict()
for (k, arr) in geo_ds.cubes
    axs = Tuple(ax isa Dim{:lon} ? geo_ds.axes[:lon] : ax for ax in arr.axes) # propagate fixed axis
    arrs[k] = YAXArray(axs, arr.data, arr.properties)
end
geo_ds = Dataset(; properties=geo_ds.properties, arrs...)
dggs2 = to_dggs_pyramid(geo_ds, level)
l2 = dggs2[2]
@test l2.area.data |> eltype == Union{Missing,Float32}
@test l2.msk_rgn.data |> eltype == Union{Missing,Int32}
@test maximum(dggs2.levels) == level
@test minimum(dggs2.levels) == 2
@test dggs2.attrs == dggs2[2].attrs
@test intersect(dggs2.attrs, dggs2[3].tas.attrs) |> length >= 0
@test length(dggs2.attrs) >= length(dggs2[2].tas.attrs)

@test dggs2[2] isa DGGSLayer
@test dggs2[level=2] isa DGGSLayer
@test dggs2[level=2, id=:tas] isa DGGSArray
@test dggs2[level=2, id=:tas, Time=1] isa DGGSArray
@test dggs2[level=2, id=:tas, Time=1].data |> size == (2, 2, 12)
@test l2[id=:tas] isa DGGSArray
@test l2[id=:tas, Time=1] isa DGGSArray
@test l2[id=:tas, Time=1].data |> size == (2, 2, 12)

@test DGGSLayer([l2.tas, l2.pr]) isa DGGSLayer
@test_throws ErrorException [p[4].tas, p[5].pr] |> DGGSLayer
@test_throws ErrorException [p[4].tas, p[4].tas] |> DGGSLayer

#
# Write pyramids 
#

d = tempname()
write_dggs_pyramid(d, dggs2)
dggs2a = open_dggs_pyramid(d)
@test dggs2.attrs == dggs2a.attrs
@test dggs2.levels == dggs2a.levels
rm(d, recursive=true)

#
# Build pyramids correctly
#

data = zeros(4, 4, 12)
data[:, :, 2] = [0 0 0 0; 0 0 0 0; 0 0 1 1; 0 0 1 1]
level = 3
axs = (Dim{:q2di_i}(1:4), Dim{:q2di_j}(1:4), Dim{:q2di_n}(1:12))
props = deepcopy(DGGS.Q2DI_DGGS_PROPS)
props["dggs_level"] = level
a = YAXArray(axs, data, props) |> DGGSArray
p = to_dggs_pyramid(a)
expected = [0 0; 0 0.625]
result = p[2].layer.data[q2di_n=2].data
@test map(ismissing, expected) == map(ismissing, result)
@test expected[2, 2] == result[2, 2]


#
# Rasters
#

a1 = Raster(rand(5, 5), (X(1:5), Y(1:5))) |> x -> to_dggs_array(x, 4; lon_name=:X, lat_name=:Y)
a2 = DimArray(rand(5, 5), (X(1:5), Y(1:5))) |> x -> to_dggs_array(x, 2; lon_name=:X, lat_name=:Y)
a3 = YAXArray((X(1:5), Y(1:5)), rand(5, 5), Dict()) |> x -> to_dggs_array(x, 2; lon_name=:X, lat_name=:Y)

@test Raster(rand(361, 181), (X(-180:180), Y(-90:90))) |> x -> to_dggs_layer(x, 4; lon_name=:X, lat_name=:Y) isa DGGSLayer
@test Raster(rand(361, 181), (X(-180:180), Y(-90:90))) |> x -> to_dggs_pyramid(x, 4; lon_name=:X, lat_name=:Y) isa DGGSPyramid

#
# plotting
#

@test plot(a; resolution=100) isa Figure
@test plot(a; type=:globe, resolution=100) isa Figure
@test plot(a; type=:native) isa Figure
@test plot(a2; type=:map) isa Figure
@test plot(a; type=:map, longitudes=-180:5:180, latitudes=-90:5:90) isa Figure


#
# Arithmetics
#

@test (a .+ 5) isa DGGSArray
@test (a * 5) isa DGGSArray
@test (a .* 5) isa DGGSArray
@test (a .> 5) isa DGGSArray

@test length(a.attrs) >= length((a * 5).attrs) # meta data invalidated after transformation
@test map(x -> x.data.data |> collect, [a * 5, 5 * a]) |> allequal

#
# getindex 
#

p_test = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/test")
l = p_test[6]
a = p_test[level=6, id=:cos]
c = Q2DI(2, 1, 14)

@test p_test isa DGGSPyramid
@test l isa DGGSLayer
@test a isa DGGSArray

@test p_test[level=6, id=:cos, center=c, radii=1:2, Time=[1, 3]] |> size == (7, 2)

@test l.cos == l[:cos] == l[id=:cos]
@test l[id=:cos, center=c, Time=[1, 3]] |> size == (2,)
@test l[id=:cos, center=c, radii=5, Time=[1, 3]] |> size == (24, 2)
@test l[id=:cos, center=c, radii=1:5, Time=[1, 3]] |> size == (61, 2)
@test l[id=:cos, Time=[1, 3]] |> size == (32, 32, 12, 2)

@test a[c] |> size == (10,)
@test a[c, 5] |> size == (24, 10)
@test a[c, 1:5] |> size == (61, 10)

@test a[c, Time=[1, 3]] |> size == (2,)
@test a[c, 5, Time=[1, 3]] |> size == (24, 2)
@test a[c, 1:5, Time=[1, 3]] |> size == (61, 2)
@test a[Time=[1, 3]] |> size == (32, 32, 12, 2)

@test a[10, 8, 8] |> size == (10,)
@test a[10, 8, 8, 1:2] |> size == (7, 10)
@test a[10, 8, 8, 1:3] |> size == (19, 10)
@test a[10, 8, 8, 2] |> size == (6, 10)
@test a[10, 8, 8, 3] |> size == (12, 10)

@test a[11.586, 50.927] == a[lon=11.586, lat=50.927]
@test a[11.586, 50.927, 1:2] |> size == (7, 10)
@test a[11.586, 50.927, 1:2, Time=[1, 3]] |> size == (7, 2)
@test a[Time=1] |> size == (32, 32, 12)
@test a[Time=1, q2di_n=3] |> size == (32, 32)

#
# Neighbors
#

p_test = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/test")
c = Q2DI(2, 1, 14)
a = p_test[6].quads
@test p_test[level=6, id=:cos, Time=[1, 2, 3], center=Q2DI(5, 1, 18), radii=1:5] |> size == (61, 3)
@test a[-52.0978195, 49.5172566, 1:5] == a[Q2DI(2, 1, 14), 1:5]
@test length(a[c, 5, :ring]) < length(a[c, 5, :disk]) < length(a[c, 5, :window])
@test p_test[6].quads[Q2DI(2, 1, 14), 1:5] |> size == (61,)
@test p_test[6].quads[Q2DI(2, 1, 14), 2] |> unique |> sort == [2, 6]
@test p_test[6].quads[Q2DI(2, 25, 1), 1:5] |> unique |> sort == [2, 11]

@test all(p_test[6].edge_disks[Q2DI(3, 28, 1), 1:5] .== 1)
@test all(p_test[6].edge_disks[Q2DI(3, 5, 1), 1:5] .== 1)
@test all(p_test[6].edge_disks[Q2DI(3, 1, 28), 1:5] .== 1)
@test all(p_test[6].edge_disks[Q2DI(3, 32, 28), 1:5] .== 1)
@test all(p_test[6].edge_disks[Q2DI(5, 18, 1), 1:5] .== 1)
@test all(p_test[6].edge_disks[Q2DI(5, 18, 32), 1:5] .== 1)
@test all(p_test[6].edge_disks[Q2DI(5, 32, 12), 1:5] .== 1)
@test all(p_test[6].edge_disks[Q2DI(5, 1, 18), 1:5] .== 1)

#
# setindex
#

lon_range = X(-180:180)
lat_range = Y(-90:90)
time_range = Ti(1:10)
level = 6
data = [exp(cosd(lon)) + t * (lat / 90) for lon in lon_range, lat in lat_range, t in time_range]
geo_arr = YAXArray((lon_range, lat_range, time_range), data, Dict())
a = to_dggs_array(geo_arr, level; lon_name=:X, lat_name=:Y)

a[Q2DI(2, 10, 10)] .= 5
a[Q2DI(3, 10, 10), Ti=1] = 5

@test all(collect(a[Q2DI(2, 10, 10)]) .== 5)
@test collect(a[Q2DI(3, 10, 10), Ti=1])[1] == 5
@test collect(a[Q2DI(4, 10, 10), Ti=1])[1] != 5
