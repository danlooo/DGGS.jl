# Convert NetCDF Climate Data

This tutorial shows how to convert climate data into a DGGS.

Download and load the [unidata example file](https://www.unidata.ucar.edu/software/netcdf/examples/files.html) from the Community Climate System Model (CCSM), one time step of precipitation flux, air temperature, and eastward wind:

```@example ncdf
using YAXArrays
using DimensionalData
using NetCDF
using Downloads

url = "https://archive.unidata.ucar.edu/software/netcdf/examples/sresa1b_ncar_ccsm3-example.nc"
path = Downloads.download(url, tempname() * ".nc")
geo_ds = open_dataset(path)
```

Plot temperature at the first time point:

```@example ncdf
using GLMakie
plot(geo_ds.tas[time=1])
```

Note, that the longitude dimension must be shifted to fit within [-180,180]:

```@example ncdf
function shift_lon(arr)
    old_lon = arr.lon.val
    n = length(old_lon) ÷ 2
    new_lon = vcat(old_lon[n+1:end] .- 360, old_lon[1:n])
    data_reordered = circshift(collect(arr.data), (n, 0, 0))
    arr2 = YAXArray((new_lon |> X, arr.lat.val |> Y, arr.time), data_reordered, arr.properties)
    return arr2
end

arr_tas = shift_lon(geo_ds.tas)
arr_pr = shift_lon(geo_ds.pr)
geo_ds2 = Dataset(tas = arr_tas, pr = arr_pr)
```

Convert to DGGS:

```@example ncdf
using DGGS
p = to_dggs_pyramid(geo_ds2, 5, "EPSG:4326")
```

Plot DGGS data at a given spatial resolution, layer and time:

```@example ncdf
plot(p[5].tas[time=1])
```
