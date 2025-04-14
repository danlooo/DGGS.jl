using DGGS
using YAXArrays
using DimensionalData
using BenchmarkTools
using Printf

@info "Start benchmark using $(nthreads()) threads"

resolution = 11
# download tif tile here: https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/GRAY_50M_SR.zip
geo_ds = open_dataset("etc/naturalearth_50m_gray_earth.tif")
geo_coords = [(lon, lat) for lon in geo_ds.X for lat in geo_ds.Y]
@info @printf "Benchmark projection of %e geo coordinates to cells (50m resolution)" length(geo_coords)
@benchmark to_cell(geo_coords, resolution)

cells = [Cell(i, j, n, resolution) for n in 0:4 for j in 0:2^resolution-1 for i in 0:2*2^resolution-1]
@info @printf "Benchmark reprojection of %e cells to geo coords (resolution %i i.e. 50m)" length(cells) resolution
@benchmark to_geo(cells)