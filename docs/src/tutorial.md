# Tutorial

## Create a Discrete Global Grid System (DGGS)

Create a data cube in geographical coordinates:

```@example dggs
using YAXArrays
lon_range = -180:180
lat_range = -90:90
data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
axlist = [
    RangeAxis("lon", lon_range),
    RangeAxis("lat", lat_range)
]
geo_cube = YAXArray(axlist, data)
```

Let's create a DGGS using Synder Equal Area projection (`"ISEA"`), an aperture of 4 (number of child cells of a given parent cell), a hexagonal grid shape at 3 different resolutions:

```@example dggs
using DGGS
dggs = GridSystem(geo_cube, "ISEA", 4, "HEXAGON", 3)
```

The data cube at the highest resolution has only one spatial index dimension, i.e. the cell id:

```@example dggs
get_cell_cube(dggs, 3)
```

Plotting requires re-griding to geographical coordinates:

```@example dggs
using CairoMakie
using GeoMakie
dggs_geo_cube = get_geo_cube(dggs, 3)
fig = Figure(resolution=(1000, 1000))
ga1 = GeoAxis(fig[1, 1]; dest="+proj=ortho", coastlines=true, lonlims=(-90, 90))
surface!(ga1, lon_range, lat_range, dggs_geo_cube.data; colormap=:rainbow_bgyrm_35_85_c69_n256, shading=false)
fig
```

## Explore the grid

Let's create our first grid from a preset:

```@example 1
using DGGS
grid = create_toy_grid()
```

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

Convert points between cell id and geographic coordinates:

```@example 1
get_cell_ids(grid, 80, -170)
```

and vice versa:

```@example 1
get_geo_coords(grid, 5)
```

## Convert a data cube

A data cube is an n-dimensional array in which we have a value for each possible combination of indices, e.g., a temperature value for each geographical coordinate and also for each time point.
Here, [YAXArrays](https://juliadatacubes.github.io/YAXArrays.jl/dev/) is used to represent the data cube.

Create a data cube with geographical coordinates using YAXArrays:

```@example 2
using YAXArrays, NetCDF
using Downloads
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
filename = Downloads.download(url, "tos_O1_2001-2002.nc") # you pick your own path
geo_cube = Cube(filename)
geo_cube
```

Indeed, we have both longitude and latitude as spatial index dimensions.
Now we can define a grid and create a new data cube `cell_cube` having just the cell id as a single spatial index dimension in accordance to the created grid:

```@example 2
using DGGS
grid = Grid("ISEA", 4, "HEXAGON", 3)
cell_cube = get_cell_cube(grid, geo_cube, "lat", "lon")
cell_cube
```

Vice versa, we can also transform a cell cube back to a geographical one:

```@example 2
geo_cube_2 = get_geo_cube(grid, cell_cube)
geo_cube_2
```

## Plot a cell data cube

The cell data cube can be visualized with the cells plotted as a choropleth map.
Cell boundaries can be drawn as shapes.
Lets start by importing a cell data cube from a NetCDF file:

```@example 3
# create the grid
using DGGS
grid = create_toy_grid()
boundaries = get_cell_boundaries(grid)

# create the cell cube
using YAXArrays, NetCDF
using Downloads
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
filename = Downloads.download(url, "tos_O1_2001-2002.nc") # you pick your own path
geo_cube = Cube(filename)
cell_cube = get_cell_cube(grid, geo_cube, "lat", "lon")
```

Now we can export the cell boundary polygons and use [GeoMakie](https://geo.makie.org/stable/) for plotting:

```@example 3
using GeoDataFrames, GeoJSON, GeoMakie, CairoMakie

boundaries_path = "boundaries.geojson"
GeoDataFrames.write(boundaries_path, boundaries)
boundaries_fc = GeoJSON.read(read(boundaries_path))

fig = Figure(resolution = (1200,800), fontsize = 22)

ax = GeoAxis(
    fig[1,1];
    dest = "+proj=wintri",
    title = "DGGS",
    tellheight = true,
)

hm2 = poly!(
    ax, boundaries_fc;
    color = cell_cube.data,
    colormap = Reverse(:plasma),
    strokecolor = :blue,
    strokewidth = 0.25
)
fig
```