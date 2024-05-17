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

p = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/sresa1b_ncar_ccsm3-example")
l = p[4]
a = l.tas

@test a.id == :tas
@test length(p.attrs) == length(l.attrs)
@test length(a.attrs) > length(p.attrs)
@test (setdiff(p.attrs, l.attrs) .|> x -> x.first) == ["_DGGS"] # same global attrs expect DGGS level


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
geo_ds.axes[:lon] = X(geo_ds.axes[:lon] .|> x -> x < 180 ? x : -(360 - x))
arrs = Dict()
for (k, arr) in geo_ds.cubes
    k == :msk_rgn && continue # exclude mask
    axs = Tuple(ax isa Dim{:lon} ? geo_ds.axes[:lon] : ax for ax in arr.axes) # propagate fixed axis
    arrs[k] = YAXArray(axs, arr.data, arr.properties)
end
geo_ds = Dataset(; properties=geo_ds.properties, arrs...)
dggs2 = to_dggs_pyramid(geo_ds, level)
@test maximum(dggs2.levels) == level
@test minimum(dggs2.levels) == 2
@test dggs2.attrs == dggs2[2].attrs
@test intersect(dggs2.attrs, dggs2[3].tas.attrs) |> length > 0
@test length(dggs2.attrs) < length(dggs2[2].tas.attrs)

#
# Write pyramids 
#

d = tempname()
write_dggs_pyramid(d, dggs2)
dggs2a = open_dggs_pyramid(d)
@test dggs2.attrs == dggs2a.attrs
@test dggs2.levels == dggs2a.levels
@test dggs2.bands == dggs2a.bands
rm(d, recursive=true)

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
@test plot(a2; type=:map) isa Figure
@test plot(a; type=:map, longitudes=-180:5:180, latitudes=-90:5:90) isa Figure