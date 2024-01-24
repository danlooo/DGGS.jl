using DGGS

lon_range = -180:180
lat_range = -90:90
level = 6
data = [exp(cosd(lon)) + 3(lat / 90) for lat in lat_range, lon in lon_range]
geo_cube = GeoCube(data, lon_range, lat_range)
cell_cube = CellCube(geo_cube, level)
dggs = GridSystem(data, lon_range, lat_range, level)

plot(dggs)

dggs_modis = GridSystem("data/modis-ndvi-2001.dggs")
query(dggs_modis[10], "Ti=2001-04-01") |> plot