
# Convert {#Convert}

## Geographical coordinates to cell ids and vice versa {#Geographical-coordinates-to-cell-ids-and-vice-versa}

Transform geographic coordinates (lon,lat) into cell ids:

```julia
using DGGS

level = 6
geo_coords = [
    (11.586, 50.927),
    (-150, 80),
    (-150.001, 80.001)
]
cell_ids = transform_points(geo_coords, level)
```


```ansi
3-element Vector{Q2DI{Int64}}:
 Q2DI(3,2,30)
 Q2DI(1,2,12)
 Q2DI(1,2,12)
```


And convert them back:

```julia
transform_points(cell_ids, level)
```


```ansi
3-element Vector{Tuple{Float64, Float64}}:
 (11.25, 51.0863507)
 (-147.2430432, 79.0526331)
 (-147.2430432, 79.0526331)
```


This will map any point within a cell to the corresponding cell center, resulting in the same output for the latter two nearby points. Manual coordinate conversion is usually not necessary, because one could use geographical coordinates directly in `getindex` methods.

## `AbstractDimArray` to `DGGSArray` {#AbstractDimArray-to-DGGSArray}

Turn any `AbstractDimArray` like a `Array`, `YAXArray`, and `Raster` into a `DGGSArray`:

```julia
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


```ansi
[37mDGGSArray{Union{Missing, Float64}, 6}[39m
Name:		layer
DGGS:		DGGRID ISEA4H Q2DI ⬢ at level 6
Attributes:	11
[37mNon spatial axes:[39m
  [31mTi[39m 10 Int64 points

```


All geographic input coordinates must be in decimal degrees in the WGS84 format. Longitudes must be within [-180,180] and latitudes must be within [-90,90]. Arrays must be in ascending order, i.e. sorted from west to east and from south to north.

This requires some adjustments in some datasets:

```julia
using NetCDF
download("https://github.com/danlooo/DGGS.jl/raw/main/test/sresa1b_ncar_ccsm3-example.nc", "example.nc")
geo_ds = open_dataset("example.nc")
geo_ds.lon, geo_ds.lat
```


```ansi
([38;5;209m↓ [39m[38;5;209mlon[39m [38;5;209m0.0f0:1.40625f0:358.59375f0[39m,
[38;5;32m→ [39m[38;5;32mlat[39m [38;5;32m-88.927734f0:1.4004368f0:88.927734f0[39m)
```


Let&#39;s change the lon axes to the desired format:

```julia
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


```ansi
[37mDGGSPyramid[39m
DGGS: DGGRID ISEA4H Q2DI ⬢
Levels: [2, 3, 4, 5, 6]
[37mNon spatial axes:[39m
  [31mtime[39m 1 CFTime.DateTimeNoLeap points
  [31mplev[39m 17 Float64 points
[37mArrays:[39m
  [31mtas[39m [37mair_temperature[39m [37m(:[39m[37mtime[39m[37m) [39m[34mK[39m Union{Missing, Float32} aggregated
  [31mua[39m [37meastward_wind[39m [37m(:plev, :time)[39m[37m [39m[34mm s-1[39m Union{Missing, Float32} aggregated
  [31mpr[39m [37mprecipitation_flux[39m [37m(:[39m[37mtime[39m[37m) [39m[34mkg m-2 s-1[39m Union{Missing, Float32} aggregated
  [31marea[39m [34mmeter2[39m Union{Missing, Float32} 

```


If an axis is in ascending order (e.g. latitude ranging from 90 to -90), it must be reversed beforehand.

## `DGGSArray` to `DGGSPyramid` {#DGGSArray-to-DGGSPyramid}

Add all coarser resolution levels to an existing `DGGSArray`:

```julia
p = to_dggs_pyramid(a)
```


```ansi
[37mDGGSPyramid[39m
DGGS: DGGRID ISEA4H Q2DI ⬢
Levels: [2, 3, 4, 5, 6]
[37mNon spatial axes:[39m
  [31mTi[39m 10 Int64 points
[37mArrays:[39m
  [31mlayer[39m [37m(:[39m[37mTi[39m[37m) [39m Union{Missing, Float64} 

```


A strided hexagonal convolution is applied to aggregate a cell and its direct neighbors to a coarser resolution having cells that are 4 times larger. The coarsest resolution is level 2 in which the surface of the earth is divided into 42 cells.

## `DGGSPyramid` to `DGGSLayer` and `DGGSArray` {#DGGSPyramid-to-DGGSLayer-and-DGGSArray}

These operations just select subsets of the pyramid and do not modify the spatial axes. See select for more details.

```julia
l = p[5]
```


```ansi
[37mDGGSLayer{5}[39m
DGGS:	DGGRID ISEA4H Q2DI ⬢ at level 5
[37mNon spatial axes:[39m
  [31mTi[39m 10 Int64 points
[37mArrays:[39m
  [31mlayer[39m [37m(:[39m[37mTi[39m[37m) [39m Union{Missing, Float64} 

```


```julia
a = p[5].layer
```


```ansi
[37mDGGSArray{Union{Missing, Float64}, 5}[39m
Name:		layer
DGGS:		DGGRID ISEA4H Q2DI ⬢ at level 5
Attributes:	11
[37mNon spatial axes:[39m
  [31mTi[39m 10 Int64 points

```

