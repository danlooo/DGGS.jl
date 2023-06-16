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

```@example dggs
using DGGS
using CairoMakie
plot_geo_cube(geo_cube)
```

Let's create a DGGS using Synder Equal Area projection (`"ISEA"`), an aperture of 4 (number of child cells of a given parent cell), a hexagonal grid shape at 3 different resolutions:

```@example dggs
dggs = GridSystem(geo_cube, "ISEA", 4, "HEXAGON", 3)
```

The data cube at the highest resolution has only one spatial index dimension, i.e. the cell id:

```@example dggs
get_cell_cube(dggs, 3)
```

Plot the DGGS at a given resolution

```@example dggs
plot_grid_system(dggs, 3)
```

A DGGS cell represent all points within its boundary polygon.
This acts as a pooling mechanism.
The hexagonal topology is easily recognizable at this low resolution.

## Explore the grid

A DGGS consists of multiple grids with varying resolutions.
Let's create our first grid to explore its properties:

```@example grid
using DGGS
grid = Grid("ISEA", 4, "HEXAGON", 2)
```

Create a data frame containing the center point or boundary polygon for all cells of the grid:

```@example grid
boundaries = get_cell_boundaries(grid)
centers = get_cell_centers(grid)

println(boundaries[1:5,:])
```

The data frames can be saved e.g. as geojson to be used in other tools like GIS:

```@example grid
using GeoDataFrames
GeoDataFrames.write("boundaries.geojson", boundaries)
```

Convert points between cell id and geographic coordinates:

```@example grid
get_cell_ids(grid, 80, -170)
```

and vice versa:

```@example grid
get_geo_coords(grid, 5)
```

The coordinates may differ slightly, because a cell covers all points of a given area and only the center point is returned.
This tutorial uses very low resolutions to demonstrate the properties of a DGGS.
In practice, much higher resolution levels should be chosen for spatial analysis, diminishing these inaccuracies.

## Import NetCDF files into a DGGS 

Here we will import data with a geographical grid into a DGGS.
Here we will explore Sea surface temperatures collected by PCMDI for use by the IPCC stored in a NetCDF file.
First, we need to create a geographical data cube.

Create a data cube with geographical coordinates using YAXArrays:

```@example netcdf
using YAXArrays
using NetCDF
using Downloads
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
filename = Downloads.download(url, "tos_O1_2001-2002.nc")
geo_cube_raw = Cube("tos_O1_2001-2002.nc")
```

The longitudes range from ~0° to ~365°.
We need to convert them to the interval -180° to 180°:

```@example netcdf
data = reverse(geo_cube_raw.data[:, :, 1]; dims = 2)
latitudes = reverse(geo_cube_raw.lat)
longitudes = geo_cube_raw.lon .- 180
axlist = [
    RangeAxis("lon", longitudes),
    RangeAxis("lat", latitudes)
]
geo_cube = YAXArray(axlist, data)
```

Plot the imported geo data cube:

```@example netcdf
using DGGS
plot_geo_cube(geo_cube)
```

Transform it into a DGGS:

```@example netcdf
dggs = GridSystem(geo_cube, "ISEA", 4, "HEXAGON", 3)
plot_grid_system(dggs, 3)
```