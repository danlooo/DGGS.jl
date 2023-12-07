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




cell_cubes = [CellCube(geo_cube, x) for x in 1:6]
cell_cubes[6].data[q2di_n=3].data |> heatmap
cell_cubes[5].data[q2di_n=3].data |> heatmap
cell_cubes[4].data[q2di_n=3].data |> heatmap
cell_cubes[3].data[q2di_n=3].data |> heatmap

# aggregate

cell_cube = cell_cubes[4]

function aggregate_cell_cube(xout, xin; agg_func=filter_null(mean))
    fac = ceil(Int, size(xin, 1) / size(xout, 1))
    for j in axes(xout, 2)
        for i in axes(xout, 1)
            iview = ((i-1)*fac+1):min(size(xin, 1), (i * fac))
            jview = ((j-1)*fac+1):min(size(xin, 2), (j * fac))
            xout[i, j] = agg_func(view(xin, iview, jview))
        end
    end
end

coarser_cell_array = mapCube(
    aggregate_cell_cube,
    cell_cube.data,
    indims=InDims(:q2di_i, :q2di_j),
    outdims=OutDims(
        Dim{:q2di_i}(range(0; step=1, length=2^(cell_cube.level - 2))),
        Dim{:q2di_j}(range(0; step=1, length=2^(cell_cube.level - 2)))
    )
)

# make pyramid
# for each quad: Aggregate 4 values to one like on a normal pyramid

cell_cube_pyramid = Vector{CellCube}(undef, cell_cube.level)
cell_cube_pyramid[cell_cube.level] = cell_cube

