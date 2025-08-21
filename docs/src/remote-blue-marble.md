# Remote data

access the remote DGGS pyramid [Blue Marble](https://github.com/danlooo/blue-marble.dggs.zarr):

```@example remote_blue_marble
using DGGS
using Zarr
using GLMakie

store = zopen("https://raw.githubusercontent.com/danlooo/blue-marble.dggs.zarr/refs/heads/master")
p = open_dggs_pyramid(store)
```

Plot the pyramid:

```@example remote_blue_marble
plot(p, :Red, :Green, :Blue; scale_factor=1/255)
```

Denkark is covered by two different faces of the icosahedron, yielding different oritentations of the DGGS zones:


```@example remote_blue_marble
using Extents
bbox_dk = Extent(X = (8.5,9.5), Y=(56,57.0))
plot(p, :Red, :Green, :Blue; scale_factor=1/255, extent=bbox_dk)
```


Subset a dataset with arrays at the same spatial resolution:

```@example remote_blue_marble
ds = p[8]
```

Extract a band:

```@example remote_blue_marble
a = ds.Blue
```