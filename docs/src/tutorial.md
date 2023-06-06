# Tutorial

## Create a grid

Let's create our first grid from a preset:

```@example 1
using DGGS
grid = create_toy_grid()
```

## Explore the grid

Create a data frame containing the center point or boundary polygon for all cells of the grid:

```@example 1
boundaries = get_cell_boundaries(grid)
centers = get_cell_centers(grid)

println(boundaries[1:5,:])
```

The data frames can be saved e.g. as geojson to be used in other tools like GIS:

```julia
using GeoDataFrames
GeoDataFrames.write("boundaries.geojson", boundaries)
```

## Convert points

Convert between cell id and geographic coordinates:

```@example 1
get_cell_ids(grid, 80, -170)
```

and vice versa:

```@example 1
get_geo_coords(grid, 5)
```

## Convert data cube

A data cube is an n-dimensional array in which we have a value for each possible combination of indices, e.g., a temperature value for each geographical coordinate and also for each time point.
Here, [YAXArrays](https://juliadatacubes.github.io/YAXArrays.jl/dev/) is used to represent the data cube.

Create a data cube with geographical coordinates using YAXArrays:

```@example 2
using YAXArrays, NetCDF
using Downloads
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
filename = Downloads.download(url, "tos_O1_2001-2002.nc") # you pick your own path
geo_cube = Cube(filename)
```

Indeed, we have both longitude and latitude as spatial index dimensions.
Now we can define a grid and create a new data cube `cell_cube` having just the cell id as a single spatial index dimension in accordance to the created grid:

```@example 2
using DGGS
grid = Grid("ISEA", 4, "HEXAGON", 3)
cell_cube = get_cell_cube(grid, geo_cube, "lat", "lon")
```