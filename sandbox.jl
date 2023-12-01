# Use YAXArray to cache cell_ids of tiles
#
# We want to cache cell_ids (Q2DI) for tile ids (xyz).
# This will be ~0.4MB per tile or ~26 GB at zoom level 8 (~600m at equator)

using Infiltrator
using BenchmarkTools

using DGGS
using YAXArrays

max_z = 3
cell_ids = Array{Matrix{DGGS.Q2DI},3}(undef, 2^max_z, 2^max_z, max_z + 1) # x,y,z

for tile_id in eachindex(DGGS.cell_ids)
    cell_ids[tile_id[1]+1, tile_id[2]+1, tile_id[3]+1] = DGGS.cell_ids[tile_id]
end

c = YAXArray(cell_ids)


@benchmark DGGS.transform_points(range(-180, 180; length=2), range(-90, 90; length=2), 6)
@benchmark DGGS.transform_points(range(-180, 180; length=15), range(-90, 90; length=15), 6)
# Overhead: time to transform 1 to ~250 geo points is const.
# A XYZ tile has 65536 points

## Explore caching
# - used for tile_id to cell_ids and to reduce the data cube to one non-spatial dimension
# - cache on disk must be much bigger than menory
# - disk cache is fasster than re-run, because we save csv write and parsing time

# Evaluation:
# - Cache can't be re-used in other sessions.
# - SimpleCaching makes a new file for each call
#
# Result:
# - Cache using a YAXArray 
using Caching
using Serialization


@cache transform_points_cached = (x, y, z, level) -> DGGS.transform_points(x, y, z, level) "data/transform_points_cached.cache.bin"
serialize("data/transform_points_cached.bin", transform_points_cached)

Threads.@threads for x in 1:10
    transform_points_cached((x, 0, 4, 6))
end

@syncache! transform_points_cached

transform_points_cached((1, 0, 4, 6))
transform_points_cached((2, 0, 4, 6))

@persist! transform_points_cached

# in another session
transform_points_cached = deserialize("data/transform_points_cached.bin", Cache)
transform_points_cached((x, 0, 4, 6))



# Load real modis data
using YAXArrays
using DimensionalData

cell_cube = CellCube("/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d050_monthly/MODIS/MOD13C2.006/Data/NDVI/NDVI.7200.3600.2001.nc", "longitude", "latitude", 6)

lon_range = -180:1.1:180
lat_range = -90:0.9:90
time_range = 1
geo_data = [t * exp(cosd(lon + (t * 10))) + 3((lat - 50) / 90) for lon in lon_range, lat in lat_range, t in time_range]
axlist = (
    Dim{:lon}(lon_range),
    Dim{:lat}(lat_range)
)
geo_array = YAXArray(axlist, geo_data)
geo_cube = GeoCube(geo_array)
cell_cube = CellCube(geo_cube, 6)
# savecube(cell_cube.data, "data/example.cube.zarr"; driver=:zarr)





using DGGS
using ProgressMeter
using ThreadSafeDicts

lons = -180:0.5:180
lats = -90:90
pixels = Iterators.product(lons, lats)
chunks = Iterators.partition(pixels, 500) |> collect
res = ThreadSafeDict()
p = Progress(length(lon_chunks))
Threads.@threads for i in eachindex(chunks)
    @infiltrate
    # DGGS.transform_points(pixels[], lats, 6)
    next!(p)
end
finish!(p)


Threads.@threads for i in 1:1e3
    run(`ls`)
end