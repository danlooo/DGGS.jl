# Tutorial

## Create a Discrete Global Grid System (DGGS)

Create a data cube in geographical coordinates:

```@example dggs
using DGGS
lon_range = -180:180
lat_range = -90:90
data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
geo_cube = GeoCube(data, lat_range, lon_range)
```

```@example dggs
using CairoMakie
plot_map(geo_cube)
```

Let's create a DGGS using Synder Equal Area projection (ISEA), an aperture of 4 (number of child cells of a given parent cell), a hexagonal grid shape at 3 different levels:

```@example dggs
dggs = DgGlobalGridSystem(geo_cube, 3, :isea, 4, :hexagon)
```

The data cube at the highest level has only one spatial index dimension, i.e. the cell id:

```@example dggs

```

Plot the DGGS at a given level

```@example dggs
plot_map(dggs[3])
```

A DGGS cell represent all points within its boundary polygon.
This acts as a pooling mechanism.
The hexagonal topology is easily recognizable at this low level.

## Create grids

Using DGGRID:

```@example grids_create
using DGGS
grid1 = DgGrid(:isea, 4, :hexagon, 3)
```

Using a vector of geographical coordinates for center points:

```@example grids_create
center_points = [-170 -80; -165.12 81.12; -160 90]
grid2 = Grid(center_points')
```

## Explore the grid

A DGGS consists of multiple grids with varying levels.
Let's create our first grid to explore its properties:

```@example grid
using DGGS
grid = DgGrid(:isea, 4, :hexagon, 3)
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
This tutorial uses very low levels to demonstrate the properties of a DGGS.
In practice, much higher level levels should be chosen for spatial analysis, diminishing these inaccuracies.


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
geo_cube_raw = YAXArrays.Cube("tos_O1_2001-2002.nc")
```

Lets have a look at the first time point of the raw data:

```@example netcdf
using CairoMakie
heatmap(geo_cube_raw[:,:,1])
```

The map is centered at the Pacific Ocean and has longitudes ranging from ~0° to ~365°.
We need to convert the map into a null meridian centric one with longitudes ranging from -180° to 180°.
One pixel represent 2° of the longitudinal axis, because there are 180 elements.
Therefore, we need to shift the data matrix by 180/2=90 pixels.

```@example netcdf
using DGGS
data = circshift(geo_cube_raw[:,:,1], 90)
latitudes = geo_cube_raw.lat
longitudes = geo_cube_raw.lon .- 180
geo_cube = GeoCube(data, latitudes, longitudes)
```

Plot the imported geo data cube:

```@example netcdf
plot_map(geo_cube)
```

Transform it into a DGGS:

```@example netcdf
dggs = DgGlobalGridSystem(geo_cube, 3, :isea, 4, :hexagon)
```

```@example netcdf
plot_map(dggs[3])
```

Since this dataset is about ocean temperature, we do not have cells on the land area.