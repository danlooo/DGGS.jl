# Get Started

## Why DGGS.jl ?

Discrete Global Grid Systems (DGGS) tessellate the surface of the earth with hierarchical cells of equal area, minimizing distortion and loading time of large geospatial datasets, which is crucial in spatial statistics and building Machine Learning models.
DGGS native data cubes use coordinate systems other than longitude and latitude to represent the special geometry of the grid needed to reduce the distortions of the individual cells or pixels.

## Installation

::: info
Currently, we only develop and test this package on Linux machines with the latest stable release of Julia.

:::

Install the latest version from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/danlooo/DGGS.jl.git")
```

## Convert coordinates

Convert into DGGS zone ids:

```@example get_started
using DGGS
lon, lat = (50.92, 11.58)
resolution = 20
cell = to_cell(lon, lat, resolution)
```

and back:

```@example get_started
to_geo(cell)
```

This will return the geographical coordinates of cell center.
The higher the resolution, the less will be the discretization error.

## Convert data into DGGS

Lets create some data in geographical space:

```@example get_started
using DimensionalData
using YAXArrays
lon_range = X(180:-1:-180)
lat_range = Y(90:-1:-90)
geo_data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
geo_array = YAXArray((lon_range, lat_range), geo_data)
```

Plot the geo data:

```@example get_started
using GLMakie
plot(geo_array)
```

Convert it into DGGS:

```@example get_started
resolution = 4
dggs_pyramid = to_dggs_pyramid(geo_array, resolution, "EPSG:4326")
```

Extract a single variable at a given spatial refinement level:

```@example get_started
dggs_array = dggs_pyramid[3].layer1
```

Plot the array:

```@example get_started
plot(dggs_array)
```

The resolution is set extremely low to demonstrate the cell shapes.
In practice, one sets it high enough to prevent losing spatial resolution.
