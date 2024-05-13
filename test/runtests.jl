using DGGS
using DimensionalData
using YAXArrays
using GLMakie
using Test
using Rasters

#
# Load and write arrays, layers, and pyramids
#

base_path = "https://s3.bgc-jena.mpg.de:9000/dggs/sresa1b_ncar_ccsm3-example"
arr = open_array("$base_path/3/tas")
l = open_layer("$base_path/3")
dggs = open_pyramid("$base_path")

@test arr.attrs |> length == 31
@test l.attrs |> length == 19
@test dggs.attrs |> length == 19
@test (setdiff(dggs.attrs, l.attrs) .|> x -> x.first) == ["_DGGS"] # same global attrs expect DGGS level

d = tempname()
write_pyramid(d, dggs)
dggs_2 = open_pyramid(d)
@test dggs.attrs == dggs_2.attrs
@test dggs.levels == dggs_2.levels
@test dggs.bands == dggs_2.bands
rm(d, recursive=true)

#
# Convert lat/lon rasters into a DGGS
#

lon_range = X(-180:180)
lat_range = Y(-90:90)
level = 6
data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
raster = DimArray(data, (lon_range, lat_range))
arr = to_array(raster, level)
@test arr.level == level

# load netcdf geo raster
# reformat lon axes from [0,360] to [-180,180]
# skip mask
geo_ds = open_dataset("data/sresa1b_ncar_ccsm3-example.nc")
geo_ds.axes[:lon] = X(geo_ds.axes[:lon] .- 180)
arrs = Dict()
for (k, arr) in geo_ds.cubes
    k == :msk_rgn && continue
    axs = Tuple(ax isa Dim{:lon} ? geo_ds.axes[:lon] : ax for ax in arr.axes)
    arrs[k] = YAXArray(axs, arr.data, arr.properties)
end
geo_ds = Dataset(; properties=geo_ds.properties, arrs...)
dggs2 = to_pyramid(geo_ds, level)
@test maximum(dggs2.levels) == level
@test minimum(dggs2.levels) == 2