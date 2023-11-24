
#
# run test
#

lon_range = -180:180
lat_range = -90:90
time_range = 1
geo_data = [t * exp(cosd(lon + (t * 10))) + 3((lat - 50) / 90) for lon in lon_range, lat in lat_range, t in time_range]
# TODO: Allow for other dimensions. But how to cretae function e.g. cell_cube[time 0 At(0)] without Meta.parse?
axlist = (
    Dim{:lon}(lon_range),
    Dim{:lat}(lat_range)
)
geo_array = YAXArray(axlist, geo_data)
geo_cube = GeoCube(geo_array)
cell_cube = CellCube(geo_cube)
geo_cube2 = GeoCube(cell_cube)

geo_cube.data.data |> heatmap
geo_cube2.data.data |> heatmap
cell_cube.data[q2di_n=At(1)].data |> heatmap
cell_cube.data[q2di_n=At(0)].data
cell_cube[Q2DI(8, 16, 0)].data
cell_cube[-180, -90].data
