
# Select {#Select}

## Select arrays {#Select-arrays}

Open a array for testing:

```julia
using DGGS
p = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/modis")
```


```ansi
[37mDGGSPyramid[39m
DGGS: DGGRID ISEA4H Q2DI ⬢
Levels: [2, 3, 4, 5, 6, 7, 8, 9, 10]
[37mNon spatial axes:[39m
  [31mtime[39m 216 Dates.DateTime points
[37mArrays:[39m
  [31mlst[39m [37m(:[39m[37mtime[39m[37m) [39m[34mK[39m Union{Missing, Float32} 
  [31mndvi[39m [37m(:[39m[37mtime[39m[37m) [39m[34mNDVI[39m Union{Missing, Float32} 

```


Select a ndvi at a given spatial resolution level:

```julia
p[5]
```


```ansi
[37mDGGSLayer{5}[39m
DGGS:	DGGRID ISEA4H Q2DI ⬢ at level 5
[37mNon spatial axes:[39m
  [31mtime[39m 216 Dates.DateTime points
[37mArrays:[39m
  [31mlst[39m [37m(:[39m[37mtime[39m[37m) [39m[34mK[39m Union{Missing, Float32} 
  [31mndvi[39m [37m(:[39m[37mtime[39m[37m) [39m[34mNDVI[39m Union{Missing, Float32} 

```


Select an array by its id:

```julia
p[5].ndvi
```


```ansi
[37mDGGSArray{Union{Missing, Float32}, 5}[39m
Name:		ndvi
Units:		NDVI
DGGS:		DGGRID ISEA4H Q2DI ⬢ at level 5
Attributes:	18
[37mNon spatial axes:[39m
  [31mtime[39m 216 Dates.DateTime points

```


Additional filtering by any non-spatial axes e.g. `Time` still results in a `DGGSArray`:

```julia
p[5].ndvi[Time=1]
```


```ansi
[37mDGGSArray{Union{Missing, Float32}, 5}[39m
Name:		NDVI
Units:		NDVI
DGGS:		DGGRID ISEA4H Q2DI ⬢ at level 5
Attributes:	18
[37mNon spatial axes:[39m

```


Further filtering will return a `YAXArray` instead:

```julia
p[5].ndvi[Time=1][q2di_n = 2]
```


```ansi
[90m┌ [39m[38;5;209m16[39m×[38;5;32m16[39m YAXArray{Union{Missing, Float32}, 2}[90m ┐[39m
[90m├────────────────────────────────────────────┴───────────── dims ┐[39m
  [38;5;209m↓ [39m[38;5;209mq2di_i[39m Sampled{Int64} [38;5;209m0:1:15[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m,
  [38;5;32m→ [39m[38;5;32mq2di_j[39m Sampled{Int64} [38;5;32m0:1:15[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m
[90m├────────────────────────────────────────────────────── metadata ┤[39m
  Dict{String, Any} with 18 entries:
  "dggs_radius"           => 6.37101e6
  "long_name"             => "monthly NDVI CMG 0.05 Deg Monthly NDVI"
  "dggs_rotation_azimuth" => 0
  "scale_factor"          => 0.0001
  "dggs_rotation_lon"     => 11.25
  "dggs_id"               => "DGGRID ISEA4H Q2DI"
  "dggs_polygon"          => "hexagon"
  "dggs_level"            => 5
  "_FillValue"            => -9999.0
  "units"                 => "NDVI"
  "name"                  => "NDVI"
  "dggs_polyhedron"       => "icosahedron"
  "missing_value"         => -9999.0
  "add_offset"            => 0.0
  "dggs_projection"       => "isea"
  "dggs_rotation_lat"     => 58.2825
  "dggs_aperture"         => 4
  "dggs_index"            => "Q2DI"
[90m├───────────────────────────────────────────────── loaded lazily ┤[39m
  data size: 1.0 KB
[90m└────────────────────────────────────────────────────────────────┘[39m
```


## Select cells and its neighbors {#Select-cells-and-its-neighbors}

Select a single cell using geographical coordinates (lon, lat):

```julia
a = p[6].ndvi
a[11.586, 50.927]
```


```ansi
[90m┌ [39m[38;5;209m216-element [39mYAXArray{Union{Missing, Float32}, 1}[90m ┐[39m
[90m├──────────────────────────────────────────────────┴───────────────────── dims ┐[39m
  [38;5;209m↓ [39m[38;5;209mtime[39m Sampled{Dates.DateTime} [38;5;209m[Dates.DateTime("2001-01-01T00:00:00"), …, Dates.DateTime("2018-12-01T00:00:00")][39m [38;5;244mForwardOrdered[39m [38;5;244mIrregular[39m [38;5;244mPoints[39m
[90m├──────────────────────────────────────────────────────────────────── metadata ┤[39m
  Dict{String, Any} with 18 entries:
  "dggs_radius"           => 6.37101e6
  "long_name"             => "monthly NDVI CMG 0.05 Deg Monthly NDVI"
  "dggs_rotation_azimuth" => 0
  "scale_factor"          => 0.0001
  "dggs_rotation_lon"     => 11.25
  "dggs_id"               => "DGGRID ISEA4H Q2DI"
  "dggs_polygon"          => "hexagon"
  "dggs_level"            => 6
  "_FillValue"            => -9999.0
  "units"                 => "NDVI"
  "name"                  => "NDVI"
  "dggs_polyhedron"       => "icosahedron"
  "missing_value"         => -9999.0
  "add_offset"            => 0.0
  "dggs_projection"       => "isea"
  "dggs_rotation_lat"     => 58.2825
  "dggs_aperture"         => 4
  "dggs_index"            => "Q2DI"
[90m├─────────────────────────────────────────────────────────────── loaded lazily ┤[39m
  data size: 864.0 bytes
[90m└──────────────────────────────────────────────────────────────────────────────┘[39m
```


Select the same cell using DGGS coordinates (n,i,j):

```julia
a[3,2,30]
```


```ansi
[90m┌ [39m[38;5;209m216-element [39mYAXArray{Union{Missing, Float32}, 1}[90m ┐[39m
[90m├──────────────────────────────────────────────────┴───────────────────── dims ┐[39m
  [38;5;209m↓ [39m[38;5;209mtime[39m Sampled{Dates.DateTime} [38;5;209m[Dates.DateTime("2001-01-01T00:00:00"), …, Dates.DateTime("2018-12-01T00:00:00")][39m [38;5;244mForwardOrdered[39m [38;5;244mIrregular[39m [38;5;244mPoints[39m
[90m├──────────────────────────────────────────────────────────────────── metadata ┤[39m
  Dict{String, Any} with 18 entries:
  "dggs_radius"           => 6.37101e6
  "long_name"             => "monthly NDVI CMG 0.05 Deg Monthly NDVI"
  "dggs_rotation_azimuth" => 0
  "scale_factor"          => 0.0001
  "dggs_rotation_lon"     => 11.25
  "dggs_id"               => "DGGRID ISEA4H Q2DI"
  "dggs_polygon"          => "hexagon"
  "dggs_level"            => 6
  "_FillValue"            => -9999.0
  "units"                 => "NDVI"
  "name"                  => "NDVI"
  "dggs_polyhedron"       => "icosahedron"
  "missing_value"         => -9999.0
  "add_offset"            => 0.0
  "dggs_projection"       => "isea"
  "dggs_rotation_lat"     => 58.2825
  "dggs_aperture"         => 4
  "dggs_index"            => "Q2DI"
[90m├─────────────────────────────────────────────────────────────── loaded lazily ┤[39m
  data size: 864.0 bytes
[90m└──────────────────────────────────────────────────────────────────────────────┘[39m
```



![](assets/rings-disks.png)
 A 3-ring and a 2-disk around a center cell

Select a 2-disk containing all neighboring cells that are at least k cells apart including a given center cell:

```julia
a[11.586, 50.927, 1:2]
```


```ansi
[90m┌ [39m[38;5;209m7[39m×[38;5;32m216[39m YAXArray{Union{Missing, Float32}, 2}[90m ┐[39m
[90m├────────────────────────────────────────────┴─────────────────────────── dims ┐[39m
  [38;5;209m↓ [39m[38;5;209mq2di_k[39m Sampled{Int64} [38;5;209m1:7[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m,
  [38;5;32m→ [39m[38;5;32mtime[39m Sampled{Dates.DateTime} [38;5;32m[Dates.DateTime("2001-01-01T00:00:00"), …, Dates.DateTime("2018-12-01T00:00:00")][39m [38;5;244mForwardOrdered[39m [38;5;244mIrregular[39m [38;5;244mPoints[39m
[90m├──────────────────────────────────────────────────────────── loaded in memory ┤[39m
  data size: 5.91 KB
[90m└──────────────────────────────────────────────────────────────────────────────┘[39m
```


This will introduce a new dimension `q2di_k` iterating over all neighbors. The ordering of cells within this dimension is deterministic but not further specified.

Select a 3-ring of cells having the same distance to the center cell:

```julia
a[11.586, 50.927, 3]
```


```ansi
[90m┌ [39m[38;5;209m12[39m×[38;5;32m216[39m YAXArray{Union{Missing, Float32}, 2}[90m ┐[39m
[90m├─────────────────────────────────────────────┴────────────────────────── dims ┐[39m
  [38;5;209m↓ [39m[38;5;209mq2di_k[39m Sampled{Int64} [38;5;209m1:12[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m,
  [38;5;32m→ [39m[38;5;32mtime[39m Sampled{Dates.DateTime} [38;5;32m[Dates.DateTime("2001-01-01T00:00:00"), …, Dates.DateTime("2018-12-01T00:00:00")][39m [38;5;244mForwardOrdered[39m [38;5;244mIrregular[39m [38;5;244mPoints[39m
[90m├──────────────────────────────────────────────────────────── loaded in memory ┤[39m
  data size: 10.12 KB
[90m└──────────────────────────────────────────────────────────────────────────────┘[39m
```


## long and short syntax {#long-and-short-syntax}

Selection on both spatial and non-spatial dimensions can be performed using keyword-based arguments on pyramids, ndvis, and arrays:

```julia
p[id=:ndvi, Time=1, level=5, lon=11.586, lat=50.927]
```


```ansi
[90m┌ [39m0-dimensional YAXArray{Union{Missing, Float32}, 0}[90m ┐[39m
[90m├────────────────────────────────────────────────────┴ metadata ┐[39m
  Dict{String, Any} with 18 entries:
  "dggs_radius"           => 6.37101e6
  "long_name"             => "monthly NDVI CMG 0.05 Deg Monthly NDVI"
  "dggs_rotation_azimuth" => 0
  "scale_factor"          => 0.0001
  "dggs_rotation_lon"     => 11.25
  "dggs_id"               => "DGGRID ISEA4H Q2DI"
  "dggs_polygon"          => "hexagon"
  "dggs_level"            => 5
  "_FillValue"            => -9999.0
  "units"                 => "NDVI"
  "name"                  => "NDVI"
  "dggs_polyhedron"       => "icosahedron"
  "missing_value"         => -9999.0
  "add_offset"            => 0.0
  "dggs_projection"       => "isea"
  "dggs_rotation_lat"     => 58.2825
  "dggs_aperture"         => 4
  "dggs_index"            => "Q2DI"
[90m├──────────────────────────────────────────────── loaded lazily ┤[39m
  data size: 4.0 bytes
[90m└───────────────────────────────────────────────────────────────┘[39m
```


which is equivalent to:

```julia
p[5].ndvi[Time=1][11.586, 50.927]
```


```ansi
[90m┌ [39m0-dimensional YAXArray{Union{Missing, Float32}, 0}[90m ┐[39m
[90m├────────────────────────────────────────────────────┴ metadata ┐[39m
  Dict{String, Any} with 18 entries:
  "dggs_radius"           => 6.37101e6
  "long_name"             => "monthly NDVI CMG 0.05 Deg Monthly NDVI"
  "dggs_rotation_azimuth" => 0
  "scale_factor"          => 0.0001
  "dggs_rotation_lon"     => 11.25
  "dggs_id"               => "DGGRID ISEA4H Q2DI"
  "dggs_polygon"          => "hexagon"
  "dggs_level"            => 5
  "_FillValue"            => -9999.0
  "units"                 => "NDVI"
  "name"                  => "NDVI"
  "dggs_polyhedron"       => "icosahedron"
  "missing_value"         => -9999.0
  "add_offset"            => 0.0
  "dggs_projection"       => "isea"
  "dggs_rotation_lat"     => 58.2825
  "dggs_aperture"         => 4
  "dggs_index"            => "Q2DI"
[90m├──────────────────────────────────────────────── loaded lazily ┤[39m
  data size: 4.0 bytes
[90m└───────────────────────────────────────────────────────────────┘[39m
```

