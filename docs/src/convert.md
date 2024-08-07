# Convert

## Geographical coordinates to cell ids and vice versa

- use transform_points
- just use lon,lat in getindex

## `AbstractDimArray` to `DGGSArray`

Turn any `AbstractDimArray` like a `YAXArray` and `Raster` into a `DGGSArray`:

```@example convert
using DGGS
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