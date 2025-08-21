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

A vertex of the icosahedron used for projection is near Gothenburg, Sweeden.
This results into different oritentations of the DGGS zones, depending on which polyhedral face they belong to:

```@example remote_blue_marble
using Extents
bbox = Extent(X = (10.5,12), Y=(57.5,59))
plot(p, :Red, :Green, :Blue; scale_factor=1/255, extent=bbox)
```


Subset a dataset with arrays at the same spatial resolution:

```@example remote_blue_marble
ds = p[8]
```

Extract a band:

```@example remote_blue_marble
a = ds.Blue
```