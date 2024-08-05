# DGGS.jl <img src="docs/src/assets/logo.drawio.svg" align="right" height="138" />

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://danlooo.github.io/DGGS.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://danlooo.github.io/DGGS.jl/dev/)
[![Build Status](https://github.com/danlooo/DGGS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/danlooo/DGGS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/danlooo/DGGS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/danlooo/DGGS.jl)

DGGS.jl is a Julia Package for scalable geospatial analysis using Discrete Global Grid Systems (DGGS), which tessellate the surface of the earth with hierarchical cells of equal area, minimizing distortion and loading time of large geospatial datasets, which is crucial in spatial statistics and building Machine Learning models.

## Important Note

This project is currently under intensive development.
The API is not considered stable yet.
There may be errors in some outputs.
We do not take any warranty for that.
Please test this package with caution.
Bug reports and feature requests are welcome.
Please create a [new issue](https://github.com/danlooo/DGGS.jl/issues/new) for this.

## Get Started

DGGS.jl currently only officially supports Julia 1.9 running on a 64bit Linux machine.
This package can be installed in Julia with the following commands:

```Julia
using Pkg
Pkg.add(url="https://github.com/danlooo/DGGS.jl.git")
```


```julia
using DGGS
p1 = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/example-ccsm3")
```
```
DGGSPyramid
DGGS: DGGRID ISEA4H Q2DI ⬢
Levels: Integer[5, 4, 6, 2, 3]
Non spatial axes:
  Time CFTime.DateTimeNoLeap
  plev Float64
Arrays: 
  tas air_temperature (:Time) K Union{Missing, Float32} aggregated
  ua eastward_wind (:plev, :Time) m s-1 Union{Missing, Float32} aggregated
  pr precipitation_flux (:Time) kg m-2 s-1 Union{Missing, Float32} aggregated
  area meter2 Union{Missing, Float32} 
```

Get temperature in Kelvin at the first time point of a hexagonal center cell and its 6 direct neighbors around a given coordinate at spatial resolution level 6:

```julia
p1[level=6, id=:tas, Time=1, lon=11.586, lat=50.927, radii=1:2] |> collect
```
```
7-element Vector{Union{Missing, Float32}}:
 282.60416f0
 282.3016f0
 282.34818f0
 281.85855f0
 281.90802f0
 282.06424f0
 282.39597f0
```

Create a DGGS based on a synthetic data in a geographical grid:

```julia
using DimensionalData
lon_range = X(-180:180)
lat_range = Y(-90:90)
level = 6
data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
raster = DimArray(data, (lon_range, lat_range))
p2 = to_dggs_pyramid(raster, level)
```
```
[ Info: Step 1/2: Transform coordinates
[ Info: Step 2/2: Re-grid the data
DGGSPyramid
DGGS: DGGRID ISEA4H Q2DI ⬢
Levels: Integer[5, 4, 6, 2, 3]
Non spatial axes:
Arrays: 
  layer  Union{Missing, Float64} 
```

Write DGGS data to disk and load them back:

```julia
write_dggs_pyramid("example.dggs", p2)
p2a = open_dggs_pyramid("example.dggs")
```

Access individual layers and arrays:

```julia
highest_layer = p1[6]
air_temperature = highest_layer.tas
temperature_map = collect(p1[level=2, id=:tas, Time=1])
```

Visualization:

```julia
using GLMakie
plot(p1[6].tas)
plot(p1[6].tas; type=:map, longitudes=-180:180, latitudes=-90:90)
```

## Development

This project is based on [DGGRID](https://github.com/sahrk/DGGRID).

## Funding

<p>
<a href = "https://earthmonitor.org/">
<img src="https://earthmonitor.org/wp-content/uploads/2022/04/european-union-155207_640-300x200.png" align="left" height="50" />
</a>

<a href = "https://earthmonitor.org/">
<img src="https://earthmonitor.org/wp-content/uploads/2022/04/OEM_Logo_Horizontal_Dark_Transparent_Background_205x38.png" align="left" height="50" />
</a>
</p>

This project has received funding from the [Open-Earth-Monitor Cyberinfrastructure](https://earthmonitor.org/) project that is part of European Union's Horizon Europe research and innovation programme under grant agreement No. [101059548](https://cordis.europa.eu/project/id/101059548).