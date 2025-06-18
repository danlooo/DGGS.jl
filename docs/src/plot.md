# Plot

```@example plot
using GLMakie
using DGGS
p = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/modis")
a = p[10].ndvi
```

Plotting is performed on a `DGGSArray` at a given spatial resolution level.
Selecting on non-spatial dimensions (e.g. Time) can be done later on in the interactive plot.
Image resolution can be adjusted using the `resolution` argument of the `plot` method.
Plotting requires to convert the DGGS space back to geographical coordinates.
The coordinate transformation is downloaded from a cache server instead of computed if the given resolution is available.

```@example plot
plot(a)
```

Plot as a map:

```@example plot
plot(a; type=:map)
```
