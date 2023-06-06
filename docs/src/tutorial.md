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

## Convert between cell id and geographic coordinates

```@example 1
get_cell_ids(grid, 80, -170)
```

and vice versa:

```@example 1
get_geo_coords(grid, 5)
```