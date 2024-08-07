# Select

## long syntax

Selection on both spatial and non-spatial dimensions can be performed using keyword-based arguments:


```@example start
using DGGS
p = open_dggs_pyramid("https://s3.bgc-jena.mpg.de:9000/dggs/datasets/example-ccsm3")
p[id=:tas, Time=1, level=5, lon=11.586, lat=50.927]
```



## short syntax

| code | output |
| --- | ---- |
| `a[lon,lat]` |  a `YAXArray` of data at one cell at the given geographical coordinate. Will keep all non-spatial axes, e.g. time.|
| `a[n,i,j]` |  a `YAXArray` of data at one cell at the given cell coordinate. Will keep all non-spatial axes, e.g. time.|