
# Get Started {#Get-Started}

## Why DGGS.jl ? {#Why-DGGS.jl-?}

Discrete Global Grid Systems (DGGS) tessellate the surface of the earth with hierarchical cells of equal area, minimizing distortion and loading time of large geospatial datasets, which is crucial in spatial statistics and building Machine Learning models. DGGS native data cubes use coordinate systems other than longitude and latitude to represent the special geometry of the grid needed to reduce the distortions of the individual cells or pixels.

## Installation {#Installation}

::: info Currently, we only develop and test this package on Linux machines with the latest stable release of Julia.

:::

Install the latest version from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/danlooo/DGGS.jl.git")
```


## Load a DGGS native data cube {#Load-a-DGGS-native-data-cube}

This dataset is based on the [example](https://www.unidata.ucar.edu/software/netcdf/examples/files.html) from the Community Climate System Model (CCSM), one time step of precipitation flux, air temperature, and eastward wind. 

```julia
using DGGS
p1 = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/example-ccsm3")
```


```ansi
[37mDGGSPyramid[39m
DGGS: DGGRID ISEA4H Q2DI ⬢
Levels: [2, 3, 4, 5, 6]
[37mNon spatial axes:[39m
  [31mTi[39m 1 CFTime.DateTimeNoLeap points
  [31mplev[39m 17 Float64 points
[37mArrays:[39m
  [31mmsk_rgn[39m [34mbool[39m Union{Missing, Bool} 
  [31mtas[39m [37mair_temperature[39m [37m(:[39m[37mTi[39m[37m) [39m[34mK[39m Union{Missing, Float32} aggregated
  [31mua[39m [37meastward_wind[39m [37m(:plev, :Ti)[39m[37m [39m[34mm s-1[39m Union{Missing, Float32} aggregated
  [31mpr[39m [37mprecipitation_flux[39m [37m(:[39m[37mTi[39m[37m) [39m[34mkg m-2 s-1[39m Union{Missing, Float32} aggregated
  [31marea[39m [34mmeter2[39m Union{Missing, Float32} 

```


This object contains multiple variables at several spatial resolutions. Let&#39;s select the air temperature at the highest resolution at the first time point:

```julia
a1 = p1[level=6, id=:tas, Time=1]
```


```ansi
[37mDGGSArray{Union{Missing, Float32}, 6}[39m
Name:		tas
Title:		model output prepared for IPCC AR4
Standard name:	air_temperature
Units:		K
DGGS:		DGGRID ISEA4H Q2DI ⬢ at level 6
Attributes:	40
[37mNon spatial axes:[39m

```


Plot the array as a globe:

```julia
using GLMakie
plot(a1)
```



![](assets/plot-get-started.png)


Get one hexagonal cell and all neighboring cells within a radius of `k` at a given geographical coordinate using `a[lon,lat,1:k]`:

```julia
a1[11.586, 50.927, 1:3]
```


```ansi
[90m┌ [39m[38;5;209m19-element [39mYAXArray{Union{Missing, Float32}, 1}[90m ┐[39m
[90m├─────────────────────────────────────────────────┴────── dims ┐[39m
  [38;5;209m↓ [39m[38;5;209mq2di_k[39m Sampled{Int64} [38;5;209m1:19[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m
[90m├──────────────────────────────────────────── loaded in memory ┤[39m
  data size: 76.0 bytes
[90m└──────────────────────────────────────────────────────────────┘[39m
```


## Create a DGGS native data cube {#Create-a-DGGS-native-data-cube}

Let&#39;s simulate traditional raster data with axes for longitude and latitude:

```julia
using DimensionalData
using YAXArrays

lon_range = X(-180:180)
lat_range = Y(-90:90)
time_range = Ti(1:10)
level = 6
data = [exp(cosd(lon)) + t * (lat / 90) for lon in lon_range, lat in lat_range, t in time_range]
geo_arr = YAXArray((lon_range, lat_range, time_range), data, Dict())
```


```ansi
[90m┌ [39m[38;5;209m361[39m×[38;5;32m181[39m×[38;5;81m10[39m YAXArray{Float64, 3}[90m ┐[39m
[90m├─────────────────────────────────┴───────────────────── dims ┐[39m
  [38;5;209m↓ [39m[38;5;209mX[39m Sampled{Int64} [38;5;209m-180:180[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m,
  [38;5;32m→ [39m[38;5;32mY[39m Sampled{Int64} [38;5;32m-90:90[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m,
  [38;5;81m↗ [39m[38;5;81mTi[39m Sampled{Int64} [38;5;81m1:10[39m [38;5;244mForwardOrdered[39m [38;5;244mRegular[39m [38;5;244mPoints[39m
[90m├─────────────────────────────────────────── loaded in memory ┤[39m
  data size: 4.99 MB
[90m└─────────────────────────────────────────────────────────────┘[39m
```


and convert it into a DGGS:

```julia
p2 = to_dggs_pyramid(geo_arr, level; lon_name=:X, lat_name=:Y)
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

