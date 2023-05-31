using DGGS

grid = toyGrid()
boundaries = cell_boundaries(grid)
centers = cell_centers(grid)

using GeoDataFrames
using GeoFormatTypes
using DataFrames
using ArchGDAL
Lon = -44 .+ rand(5);
Lat = -23 .+ rand(5);
df = DataFrame(geometry=ArchGDAL.createpoint.(Lon, Lat))

GDF.write("test.geojson", df, crs=GeoFormatTypes.EPSG(4326))


# export centers
df = DataFrame(geometry=ArchGDAL.createpoint.(Lon, Lat))


file = "out.geojson"

geometry = Vector{ArchGDAL.IGeometry}(undef, length(grid))
for i in eachindex(grid.data.data)
    geometry[i] = ArchGDAL.createpoint(grid.data.data[i][1], grid.data.data[i][2])
end
df = DataFrame(geometry=geometry)
GeoDataFrames.write(file, df, crs=GeoFormatTypes.EPSG(4326))



function export_centers(grid::Grid; filepath::String="centers.geojson")
    centers = cell_centers(grid)

    ArchGDAL.create(
        "centers.geojson",
        driver=ArchGDAL.getdriver("geojson")
    ) do ds
        ArchGDAL.createlayer(geom=ArchGDAL.wkbPoint) do layer
            for i in eachindex(grid.data.data)
                ArchGDAL.createfeature(layer) do f
                    ArchGDAL.setgeom!(f, ArchGDAL.createpoint(centers[i][1], centers[i][2]))
                end
            end
            ArchGDAL.copy(layer, dataset=ds)
        end
    end

    for i in eachindex(grid.data.data)
        geometry[i] = ArchGDAL.createpoint(grid.data.data[i][1], grid.data.data[i][2])
    end
    df = DataFrame(geometry=geometry)
    GeoDataFrames.write(filepath, df, crs=GeoFormatTypes.EPSG(4326))
end

grid = toyGrid()
export_centers(grid; filepath="centers.geojson")


using ArchGDAL
Lon = -44 .+ rand(5);
Lat = -23 .+ rand(5);
ArchGDAL.create(
    "centers.geojson",
    driver=ArchGDAL.getdriver("geojson")
) do ds
    ArchGDAL.createlayer(geom=ArchGDAL.wkbPoint) do layer
        for (lon, lat) in zip(Lon, Lat)
            ArchGDAL.createfeature(layer) do f
                ArchGDAL.setgeom!(f, ArchGDAL.createpoint(lat, lon))
            end
        end
        ArchGDAL.copy(layer, dataset=ds)
    end
end;