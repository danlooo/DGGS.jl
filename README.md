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

Create a data cube based on a geographical grid:

```julia
using DGGS
lon_range = -180:180
lat_range = -90:90
data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
geo_cube = GeoCube(data, lat_range, lon_range)
```
```
DGGS GeoCube
Element type: Float64
Latitude:     RangeAxis with 181 elements from -90 to 90
Longituide:   RangeAxis with 361 elements from -180 to 180
Size:         510.48 KB
```

Create a DGGS from it:

```julia
dggs = DgGlobalGridSystem(geo_cube, 3, :isea, 4, :hexagon)
```
```
DGGS DgGlobalGridSystem
Cells:   3 levels with up to 162 cells of type Float64
Grid:    DgGrid with hexagon topology, isea projection, and aperture of 4
Size:    1.69 KB
```

Checkout the [tutorial](https://danlooo.github.io/DGGS.jl/dev/tutorial/) for further examples.

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

This project has received funding from the from [Open-Earth-Monitor Cyberinfrastructure](https://earthmonitor.org/) project which is part of the European Union's Horizon Europe research and innovation programme under grant agreement No. [101059548](https://cordis.europa.eu/project/id/101059548).