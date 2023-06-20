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

Create a Discrete Global Grid System (DGGS) based on data stored in a NetCDF file:

```julia
using DGGS
using YAXArrays, NetCDF, Downloads
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
filename = Downloads.download(url, "tos_O1_2001-2002.nc")
geo_cube = Cube(filename)

dggs = GridSystem(geo_cube, "ISEA", 4, "HEXAGON", 3)
```
```
Discrete Global Grid System
Grid:   HEXAGON topology, ISEA projection, aperture of 4
Cells:  3 levels with up to 162 cells
Data:   YAXArray of type Vector{Float32} with 864 bytes
```

Checkout the [tutorial](https://danlooo.github.io/DGGS.jl/dev/tutorial/) for further examples.

## Development

This project is based on [DGGRID](https://github.com/sahrk/DGGRID).