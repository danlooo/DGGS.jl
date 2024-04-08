# Frequently Asked Questions (FAQ)

## Reverse geographical dimension

- Sometimes, the latitude can be arranged in ascending or descending order
- DGGS requires to have spatial axes in forward order


```julia
using YAXArrays, NetCDF
path = download("https://s3.bgc-jena.mpg.de:9000/dggs/test/inversed-lat.nc", "inversed-lat.nc")
c = Cube(path)
```

We need to invert the ordwer of the latitude dimensions:

```julia
c = @view c[axes(c, 1), reverse(axes(c, 2)), axes(c, 3)] 
```

Transform it into a `GeoCube` that can be further converted into a `CellCube` or `GridSystem`:

```julia
using DGGS
c = renameaxis!(c, "longitude" => :lon)
c = renameaxis!(c, "latitude" => :lat)
dggs = GeoCube(c)
```