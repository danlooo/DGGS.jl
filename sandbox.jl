using Infiltrator
using DGGS
using DimensionalData
using YAXArrays
using Plots

# MODIS land data
# @time cell_cube = CellCube("/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d050_monthly/MODIS/MOD13C2.006/Data/NDVI/NDVI.7200.3600.2001.nc", "longitude", "latitude", 6)


# Synthetic data

lon_range = -180:0.5:180
lat_range = -90:0.5:90
time_range = 1
geo_data = [t * exp(cosd(lon + (t * 10))) + 3((lat - 50) / 90) for lon in lon_range, lat in lat_range, t in time_range]
axlist = (
    Dim{:lon}(lon_range),
    Dim{:lat}(lat_range)
)
geo_array = YAXArray(axlist, geo_data)
geo_cube = GeoCube(geo_array)
cell_cube = CellCube(geo_cube, 6)
dggs = GridSystem(cell_cube)
saveGridSystem(dggs, "data/example.dggs.zarr")
dggs2 = GridSystem("data/example.dggs.zarr")



