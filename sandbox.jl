using Infiltrator
using DGGS
using DimensionalData
using YAXArrays
using Plots

# MODIS land data
@time cell_cube = CellCube("/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d050_monthly/MODIS/MOD13C2.006/Data/NDVI/NDVI.7200.3600.2001.nc", "longitude", "latitude", 12)
dggs = GridSystem(cell_cube)
saveGridSystem(dggs, "data/modis-ndvi.dggs.zarr")

# Synthetic data
lon_range = -180:0.5:180
lat_range = -90:0.5:90
time_range = 1:10
geo_data = [t * exp(cosd(lon + (t * 10))) + 3((lat - 50) / 90) for lon in lon_range, lat in lat_range, t in time_range]
axlist = (
    Dim{:lon}(lon_range),
    Dim{:lat}(lat_range),
    Dim{:time}(time_range)
)
geo_array = YAXArray(axlist, geo_data)
geo_cube = GeoCube(geo_array)
cell_cube = CellCube(geo_cube, 8)

dggs = GridSystem(cell_cube)
saveGridSystem(dggs, "data/example.dggs.zarr")
dggs2 = GridSystem("data/example.dggs.zarr")


#
using Zarr
g = zgroup(DirectoryStore("data/example2.dggs.zarr"))
for (level, cell_cube) in dggs.data
    z_arr = zcreate(eltype(cell_cube.data), g, "$level", size(dggs[level].data.data)...; compressor=Zarr.ZlibCompressor())
    z_arr = dggs[level].data.data
end

Cube("data/example2.dggs.zarr/2")

# combine data base
cache = DGGS.cache_xyz_to_q2di("data/xyz_to_q2di/"; max_z=5)
using Serialization
# DGGS.cache_xyz_to_q2di()
db = Dict()
map(readdir("data/xyz_to_q2di")) do file
    cell_ids = deserialize("data/xyz_to_q2di/$file")
    x = parse(Int, split(file, ".")[1])
    y = parse(Int, split(file, ".")[2])
    z = parse(Int, split(file, ".")[3])
    db[x, y, z] = cell_ids
end
serialize("data/xyz_to_q2di.cache.bin", db)


# Query
# "{host}/{path}/{query}/tile/{z}/{x}/{y}/tile.png"

dggs = GridSystem("data/modis-ndvi.dggs")
query(dggs[2], "all")
query(dggs[2], "Time=2001-01-01")