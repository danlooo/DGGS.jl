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

p = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/example-ccsm3.zarr")
l = p[4]
a = l.tas

@test a.id == :tas
@test length(p.attrs) == length(l.attrs)
@test length(a.attrs) > length(p.attrs)
@test (setdiff(p.attrs, l.attrs) .|> x -> x.first) == ["_DGGS"] # same global attrs expect DGGS level

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
a2 = to_dggs_array(raster, level)
raster2 = to_geo_array(a2, lon_range.val, lat_range.val)
@test a2.level == level

# low loss after converting back
@test all(abs.(Matrix(raster) .- Matrix(raster2)) .< 1.5)

# load netcdf geo raster
# reformat lon axes from [0,360] to [-180,180]
# skip mask
geo_ds = open_dataset("sresa1b_ncar_ccsm3-example.nc")
geo_ds.axes[:lon] = vcat(range(0, 180; length=128), range(-180, 0; length=128)) |> X
arrs = Dict()
for (k, arr) in geo_ds.cubes
    k == :msk_rgn && continue # exclude mask
    axs = Tuple(ax isa Dim{:lon} ? geo_ds.axes[:lon] : ax for ax in arr.axes) # propagate fixed axis
    arrs[k] = YAXArray(axs, arr.data, arr.properties)
end
geo_ds = Dataset(; properties=geo_ds.properties, arrs...)
dggs2 = to_dggs_pyramid(geo_ds, level)
l2 = dggs2[2]
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
props = Dict("_DGGS" => deepcopy(DGGS.Q2DI_DGGS_PROPS))
props["_DGGS"]["level"] = level
a = YAXArray(axs, data, props) |> DGGSArray
p = to_dggs_pyramid(a)
expected = [0 0; 0 0.625]
result = p[2].layer.data[q2di_n=2].data
@test map(ismissing, expected) == map(ismissing, result)
@test expected[2, 2] == result[2, 2]


#
# Rasters
#

a1 = Raster(rand(5, 5), (X(1:5), Y(1:5))) |> x -> to_dggs_array(x, 4)
a2 = DimArray(rand(5, 5), (X(1:5), Y(1:5))) |> x -> to_dggs_array(x, 2)
a3 = YAXArray((X(1:5), Y(1:5)), rand(5, 5), Dict()) |> x -> to_dggs_array(x, 2)

@test Raster(rand(361, 181), (X(-180:180), Y(-90:90))) |> x -> to_dggs_layer(x, 4) isa DGGSLayer
@test Raster(rand(361, 181), (X(-180:180), Y(-90:90))) |> x -> to_dggs_pyramid(x, 4) isa DGGSPyramid

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
# Neighbors
#

p_test = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/test.zarr")
c = Q2DI(2, 1, 14)
a = p_test[6].quads
@test a[-52.0978195, 49.5172566, 5] == a[Q2DI(2, 1, 14), 5]
@test length(a[c, 5, :ring]) < length(a[c, 5, :disk]) < length(a[c, 5, :window])
@test p_test[6].quads[Q2DI(2, 1, 14), 5] |> size == (61,)
@test p_test[6].quads[Q2DI(2, 1, 14), 5] |> unique == [2, 6]
@test p_test[6].quads[Q2DI(2, 25, 1), 5] |> unique == [11, 2]