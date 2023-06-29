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
dggs
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

## Import Earth System Data Cube (ESDC)

[Earth System Data Lab](https://www.earthsystemdatalab.net/) provides global data about 42 features about temperature, precipitation, among others.
The data is stored online in the [Earth System Data Cube (ESDC)](https://deepesdl.readthedocs.io/en/latest/datasets/ESDC/) and can be imported with the following commands:

```@example esdc
using DGGS
using EarthDataLab
esdc_cube = esdc(res="low")
geo_cube = GeoCube(esdc_cube)
```
Let's create a subset containing vegetation data of Asia in the year 2020:

```@example esdc
sub_geo_cube =  geo_cube[region="Asia", time=2020:2020, Variable="ndvi"]
```

Note that one can only plot a cube if the only remaining axes are longitude and latitude.
This we need to subset it to plot the first time point of the imported raster data:

```@example esdc
using CairoMakie
plot_map(sub_geo_cube[time = DateTime("2020-01-05T00:00:00")])
```

Now we can create a DGGS from it:

```@example esdc
dggs = DgGlobalGridSystem(sub_geo_cube, :isea4h, 6)
```

And plot the DGGS at the first time point with different resolutions:

```@example esdc
plot_map(dggs[5][time=DateTime("2020-01-01T01:00:00")])
```

```@example esdc
plot_map(dggs[4][time=DateTime("2020-01-01T01:00:00")])
```

Due to the nature of the vegetation index, areas in the oceans are undefined.
Thus, parent cells having any undefined value are undefined as well, resulting in undefined coastal areas at lower resolution.
We can prevent this by providing a custom aggregation function:

```@example esdc
using Statistics
aggregate_skip_undef(value) = value |> filter(! isnan) |> mean

dggs2 = DgGlobalGridSystem(sub_geo_cube, :isea4h, 6; aggregation_function = aggregate_skip_undef)
plot_map(dggs2[4][time=DateTime("2020-01-01T01:00:00")])
```

## Import Zarr Arrays

Load the zarr data into a YAXArray:

```@example zarr
using Zarr, YAXArrays
url = "gs://cmip6/CMIP6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/3hr/tas/gn/v20190710/"
geo_dataset = open_dataset(zopen(url, consolidated=true))
geo_array = geo_dataset["tas"]
```

Transform the cooordinates and create the GeoCube:

```@example zarr
using DGGS
data = circshift(geo_array[:, :, 1], length(geo_array.lon) / 2)
latitudes = geo_array.lat
longitudes = geo_array.lon .- 180
geo_cube = GeoCube(data, latitudes, longitudes)
plot_map(geo_cube)
```

## Import NetCDF files

Here we will explore Sea surface temperatures collected by PCMDI for use by the IPCC stored in a NetCDF file.
Download the NetCDF file into a YAXArray:

```@example netcdf
using YAXArrays
using NetCDF
using Downloads
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
filename = Downloads.download(url, "tos_O1_2001-2002.nc")
geo_array = YAXArrays.Cube("tos_O1_2001-2002.nc")
```

Lets have a look at the first time point of the raw data:

```@example netcdf
using CairoMakie
heatmap(geo_array[:,:,1])
```

Transform the cooordinates and create the GeoCube:

```@example netcdf
using DGGS
data = circshift(geo_array[:,:,1], 90)
latitudes = geo_array.lat
longitudes = geo_array.lon .- 180
geo_cube = GeoCube(data, latitudes, longitudes)
```

Plot the imported geo data cube:

```@example netcdf
plot_map(geo_cube)
```
Since this dataset is about ocean temperature, we do not have cells on the land area.