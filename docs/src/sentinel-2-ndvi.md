# Calculate Sentinel-2 NDVI

Load and subset geographical data:

```@example s2_ndvi
using DGGS
using YAXArrays
using ArchGDAL

geo_ds = open_dataset("https://github.com/danlooo/DGGS.jl/raw/refs/heads/main/test/data/s2-ndvi-irrigation.tif")
plot(geo_ds.nir)
```

Convert it into a DGGS array:

```@example s2_ndvi
dggs_ds = to_dggs_dataset(geo_ds, 19, geo_ds.nir.properties["projection"])
```

Calculate the NDVI:

```@example s2_ndvi
uint16_max = 65535
nir = dggs_ds.nir ./ uint16_max
red = dggs_ds.red ./ uint16_max
ndvi = @. (nir - red) / (nir + red)
```

Hereby, we use the dot assignment macro `@.` to broadcast the function over all cells individually.

Plot the NDVI in DGGS:

```@example s2_ndvi
plot(ndvi)
```
