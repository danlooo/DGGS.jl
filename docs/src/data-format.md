# Data Format

DGGS.jl stores data of individual DGGS cells in rasters, in which the cell ids can be calculated by row and column numbers.

## Projection

Geographical coordinates are first projected into a plane using Snyder equal-area projection, followed by a linear transformation, yielding a tesselation of the earth surface into 5 matrices.
The transformed coordinates are then adjusted downward to get the final discrete cell ids as (i,j,n) tuples at row 0<=i<=2*2^r, column  0<=i<=2^r, and matrix 0<=n<=4>.
Moreover, the cell id can be re-interpreted as a single integer, by concatenating the bits of the individual axial coordinates.
The refinement level r determines the spatial resolution, halving the width and length of the 5 matrices in each subsequent coarser levels.
The final DGGS in an image pyramid of datasets at a given refinement level and all coarser levels.

![](https://github.com/danlooo/DGGS.jl/raw/main/docs/src/assets/pentacube-overview.png)

## Data Format

The data is stored according to the [Unidata's Common Data Model (CDM)](https://docs.unidata.ucar.edu/netcdf-c/current/netcdf_data_model.html) using [CF Metadata Conventions](https://cfconventions.org/) by extending attributes and dimensions with the prefix `dggs_`.
An example in Zarr format is provided at [https://github.com/danlooo/blue-marble.dggs.zarr](https://github.com/danlooo/blue-marble.dggs.zarr).
The root group contains a group for all spatial refinement levels (dataset) that contain all variables (arrays) for variables and dimensions at the given level.
Each dataset MUST contain spatial dimensions `dggs_i`, `dggs_j`, and `dggs_n`.

## Meta Data

Attributes used to descirbe the DGGS:

| key | format | example |description |
| --- | --- | --- | --- |
| dggs_resolution | int | 2 | Spatial refinement level defined by the given DGGSRS |
| dggs_bbox | `{"X":[float,float],"Y":[float,float]}` | `{"X":[-180.0,180.0],"Y":[-90.0,90.0]}` | Bounding box of the data in WGS84 |
| dggs_dggsrs | string | "ISEA4D.Penta" | Name of DGGSRS |

The root group `/` MUST contain attributes `dggs_bbox` and `dggs_dggsrs`.
Dataset groups `/dggs_s{r}` MUST contain attributes `dggs_dggsrs` and `dggs_resolution`.
Arrays `/dggs_s{r}/{array_name}` MUST contain attributes `dggs_bbox`, `dggs_dggsrs`, and `dggs_resolution`.