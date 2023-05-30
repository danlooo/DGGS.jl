```@meta
CurrentModule = DGGS
```

# DGGS

Documentation for [DGGS](https://github.com/danlooo/DGGS.jl).

```@index
```

```@autodocs
Modules = [DGGS]
```

## Get Started

Install DGGS.jl:

```{julia}
using Pkg
Pkg.add("DGGS")
```

Let's create our first grid from a preset:

```{julia}
grid = toyGrid()
grid.spec

boundaries = cell_boundaries(grid)
centers = cell_centers(grid)
```