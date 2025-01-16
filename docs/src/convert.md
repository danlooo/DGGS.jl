# Convert

## Geographical coordinates to cell ids and vice versa

Transform geographic coordinates (lon,lat) into cell ids:

```@example convert
using DGGS

level = 6
geo_coords = [
    (11.586, 50.927),
    (-150, 80),
    (-150.001, 80.001)
]
cell_ids = transform_points(geo_coords, level)
```

And convert them back:

```@example convert
transform_points(cell_ids, level)
```

This will map any point within a cell to the corresponding cell center, resulting in the same output for the latter two nearby points.
Manual coordinate conversion is usually not necessary, because one could use geographical coordinates directly in `getindex` methods.

## `AbstractDimArray` to `DGGSArray`

Turn any `AbstractDimArray` like a `Array`, `YAXArray`, and `Raster` into a `DGGSArray`:

```@example convert
using DimensionalData
using YAXArrays

lon_range = X(-180:180)
lat_range = Y(-90:90)
time_range = Ti(1:10)
level = 6
data = [exp(cosd(lon)) + t * (lat / 90) for lon in lon_range, lat in lat_range, t in time_range]
geo_arr = YAXArray((lon_range, lat_range, time_range), data, Dict())
a = to_dggs_array(geo_arr, level; lon_name=:X, lat_name=:Y)
```

All geographic input coordinates must be in decimal degrees in the WGS84 format.
Longitudes must be within [-180,180] and latitudes must be within [-90,90].
Arrays must be in ascending order, i.e. sorted from west to east and from south to north.

This requires some adjustments in some datasets:

```@example convert
using NetCDF
download("https://github.com/danlooo/DGGS.jl/raw/main/test/sresa1b_ncar_ccsm3-example.nc", "example.nc")
geo_ds = open_dataset("example.nc")
geo_ds.lon, geo_ds.lat
```

Let's change the lon axes to the desired format:

```@example convert
geo_ds.axes[:lon] = vcat(range(0, 180; length=128), range(-180, 0; length=128)) |> lon
arrs = Dict()
for (k, arr) in geo_ds.cubes
    k == :msk_rgn && continue # exclude mask
    axs = Tuple(ax isa lon ? geo_ds.axes[:lon] : ax for ax in arr.axes) # propagate fixed axis
    arrs[k] = YAXArray(axs, arr.data, arr.properties)
end
geo_ds = Dataset(; properties=geo_ds.properties, arrs...)
p = to_dggs_pyramid(geo_ds, level)
```

If an axis is in ascending order (e.g. latitude ranging from 90 to -90), it must be reversed beforehand.

## `DGGSArray` to `DGGSPyramid`

Add all coarser resolution levels to an existing `DGGSArray`:

```@example convert
p = to_dggs_pyramid(a)
```

A strided hexagonal convolution is applied to aggregate a cell and its direct neighbors to a coarser resolution having cells that are 4 times larger.
The coarsest resolution is level 2 in which the surface of the earth is divided into 42 cells.

## `DGGSPyramid` to `DGGSLayer` and `DGGSArray`

These operations just select subsets of the pyramid and do not modify the spatial axes.
See select for more details.

```@example convert
l = p[5]
```

```@example convert
a = p[5].layer
```
