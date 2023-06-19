import YAXArrays: YAXArray, Cube, Cubes.formatbytes, Cubes.cubesize
import Statistics: mean
import Makie
import GeoMakie

abstract type DGGSCube end

struct CellCube <: DGGSCube
    data::YAXArray
    grid::AbstractGrid
    cell_ids

    function CellCube(data::YAXArray, grid::AbstractGrid)
        hasproperty(data, :cell_id) ? true : @error "CellCube must have property cell_id"
        eltype(data.cell_id) <: Int ? true : @error "Field cell_id must be an Integer"

        new(data, grid, data.cell_id)
    end
end

struct GeoCube <: DGGSCube
    data::YAXArray
    longitudes
    latitudes

    function GeoCube(data)
        hasproperty(data, :lon) ? true : @error "GeoCube must have property lon"
        hasproperty(data, :lat) ? true : @error "GeoCube must have property lat"
        issorted(data.lon) ? true : @error "Longitude must be sorted"
        issorted(data.lat) ? true : @error "Latitude must be sorted"
        first(data.lon) == -180 && last(data.lon) == 180 ? true : @warn "Longitude grid is not global and does not range from -180 to 180"
        first(data.lat) == -90 && last(data.lat) == 90 ? true : @warn "Latitude grid is not global and does not range from -90 to 90"

        new(data, data.lon, data.lat)
    end
end


function Base.show(io::IO, ::MIME"text/plain", geo_cube::GeoCube)
    println(io, "DGGS GeoCube")
    println(io, "Element type: $(eltype(geo_cube))")
    println(io, "Latitude:     RangeAxis with $(length(geo_cube.latitudes)) elements from $(first(geo_cube.latitudes)) to $(last(geo_cube.latitudes))")
    println(io, "Longituide:   RangeAxis with $(length(geo_cube.longitudes)) elements from $(first(geo_cube.longitudes)) to $(last(geo_cube.longitudes))")
    println(io, "size:         $(formatbytes(cubesize(geo_cube.data)))")
end

Base.eltype(geo_cube::GeoCube) = eltype(geo_cube.data)

function GeoCube(array::YAXArray, latitude_name, longitude_name)
    latitude_symbol = Symbol(latitude_name)
    longitude_symbol = Symbol(longitude_name)

    latitude_symbol in propertynames(array) ? true : error("Missing dimension $(latitude_name)")
    longitude_symbol in propertynames(array) ? true : error("Missing dimension $(longitude_name)")
    latitude_symbol != longitude_symbol ? true : error("Dimensions of longitude and latitude must be different")

    renameaxis!(array, latitude_symbol => :lat)
    renameaxis!(array, longitude_symbol => :lon)

    return GeoCube(array)
end

function GeoCube(filepath::String, latitude_name, longitude_name)
    array = YAXArrays.Cube(filepath)
    GeoCube(array, latitude_name, longitude_name)
end


"""
Export cell data cube into a traditional geographical one

Transforms a data cube with one spatial index dimensions, i. e., the cell id,
into a traditional geographical data cube with two spatial index dimensions longitude and latitude.
Values are taken from the nearest cell.
"""
function GeoCube(cell_cube::CellCube)
    longitudes = -180:180
    latitudes = -90:90

    regridded_matrix = Matrix{eltype(cell_cube)}(undef, length(longitudes), length(latitudes))

    for (lon_i, lon) in enumerate(longitudes)
        for (lat_i, lat) in enumerate(latitudes)
            cur_cell_id = get_cell_ids(cell_cube.grid, lat, lon)
            regridded_matrix[lon_i, lat_i] = cell_cube.data[cur_cell_id]
        end
    end

    axlist = [
        RangeAxis("lon", longitudes),
        RangeAxis("lat", latitudes)
    ]

    cube = YAXArray(axlist, regridded_matrix) |> GeoCube
    return cube
end

function plot!(geo_cube::GeoCube)
    # Can not use Makie plot recipies, because we need to specify the axis for GeoMakie
    # see https://discourse.julialang.org/t/accessing-axis-in-makie-plot-recipes/66006

    fig = Figure()
    ga1 = GeoAxis(fig[1, 1]; dest="+proj=wintri", coastlines=true)
    sf = surface!(ga1, geo_cube.longitudes, geo_cube.latitudes, geo_cube.data.data; colormap=:viridis, shading=false)
    cb1 = Colorbar(fig[1, 2], sf; label="Value", height=Relative(0.5))
    fig
end

"""
Import geographical data cube into a DGGS

Transforms a data cube with spatial index dimensions longitude and latitude
into a data cube with the cell id as a single spatial index dimension.
Re-gridding is done using the average value of all geographical coordinates belonging to a particular cell defined by the grid specification `grid_spec`.
"""
function CellCube(geo_cube::GeoCube, grid::AbstractGrid; aggregate_function::Function=mean)
    cell_ids = get_cell_ids(grid, geo_cube.latitudes, geo_cube.longitudes)
    cell_values = Vector{eltype(geo_cube)}(undef, length(grid))

    for cell_id in unique(cell_ids)
        cell_coords = findall(isequal(cell_id), cell_ids)
        if isempty(cell_coords)
            continue
        end
        cell_values[cell_id] = aggregate_function(geo_cube.data.data'[cell_coords])
    end

    axlist = [RangeAxis("cell_id", range(1, length(grid)))]
    cell_cube = YAXArray(axlist, cell_values)
    return CellCube(cell_cube, grid)
end

function Base.show(io::IO, ::MIME"text/plain", cell_cube::CellCube)
    println(io, "DGGS CellCube")
    println(io, "Element type: $(eltype(cell_cube))")
    println(io, "Cell id:      RangeAxis with $(length(cell_cube.cell_ids)) elements from $(first(cell_cube.cell_ids)) to $(last(cell_cube.cell_ids))")
    println(io, "size:         $(YAXArrays.Cubes.formatbytes(YAXArrays.Cubes.cubesize(cell_cube.data)))")
end

function plot!(cell_cube::CellCube)
    geo_cube = GeoCube(cell_cube)
    plot(geo_cube)
end

# Would it be better to make CellCube <: YAXArray ?
# pro: we do not need to forward implementations like this
# Cons: No abstraction layer, difficult to switch backends in the future
Base.getindex(cell_cube::CellCube, i...) = cell_cube.data[i...]
Base.firstindex(cell_cube::CellCube) = first(cell_cube.cell_ids)
Base.lastindex(cell_cube::CellCube) = last(cell_cube.cell_ids)
Base.length(cell_cube::CellCube) = length(cell_cube.cell_ids)
Base.eltype(cell_cube::CellCube) = eltype(cell_cube.data)