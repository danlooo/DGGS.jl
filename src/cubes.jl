using YAXArrays
import YAXArrays: Cubes.formatbytes, Cubes.cubesize, Cubes.getattributes, Cubes.subsetcube
import Statistics: mean
using Makie
using GeoMakie

"Subtypes must have field data::YAXArray"
abstract type AbstractCube end
cubesize(cube::AbstractCube) = YAXArrays.Cubes.cubesize(cube.data)
Base.size(cube::AbstractCube) = YAXArrays.Cubes.size(cube.data)
Base.eltype(cube::AbstractCube) = eltype(cube.data)

struct CellCube <: AbstractCube
    data::YAXArray
    grid::AbstractGrid
    cell_ids

    function CellCube(data::YAXArray, grid::AbstractGrid)
        hasproperty(data, :cell_id) || throw(ArgumentError("CellCube must have property cell_id"))
        eltype(data.cell_id) <: Int || throw(ArgumentError("Field cell_id must be an Integer"))

        new(data, grid, data.cell_id)
    end
end

struct GeoCube <: AbstractCube
    data::YAXArray
    longitudes
    latitudes

    function GeoCube(data)
        hasproperty(data, :lon) || throw(ArgumentError("GeoCube must have property lon"))
        hasproperty(data, :lat) || throw(ArgumentError("GeoCube must have property lat"))
        first(data.lon) == -180 && last(data.lon) == 180 ? true : @warn "Longitude grid is not global and does not range from -180 to 180"
        first(data.lat) == -90 && last(data.lat) == 90 ? true : @warn "Latitude grid is not global and does not range from -90 to 90"

        new(data, data.lon, data.lat)
    end
end

function Base.show(io::IO, ::MIME"text/plain", cube::AbstractCube)
    println(io, "DGGS $(typeof(cube))")
    println(io, "Element type:       $(eltype(cube))")
    println(io, "Size:               $(formatbytes(cubesize(cube)))")
    println(io, "Axes:")
    for axis in cube.data.axes
        println(io, repr(axis))
    end
    foreach(YAXArrays.Cubes.getattributes(cube.data)) do p
        if p[1] in ("labels", "name", "units")
            println(io, p[1], ": ", p[2])
        end
    end
end

function Base.getindex(geo_cube::GeoCube; kwargs...)
    subset_geo_cube = YAXArrays.Cubes.subsetcube(geo_cube.data; kwargs...)
    GeoCube(subset_geo_cube)
end

function GeoCube(array::YAXArray, latitude_name, longitude_name)
    latitude_symbol = Symbol(latitude_name)
    longitude_symbol = Symbol(longitude_name)

    latitude_symbol in propertynames(array) || throw(ArgumentError("Missing dimension $(latitude_name)"))
    longitude_symbol in propertynames(array) || throw(ArgumentError("Missing dimension $(longitude_name)"))
    latitude_symbol != longitude_symbol || throw(ArgumentError("Dimensions of longitude and latitude must be different"))

    renameaxis!(array, latitude_symbol => :lat)
    renameaxis!(array, longitude_symbol => :lon)

    return GeoCube(array)
end

function GeoCube(filepath::String, latitude_name, longitude_name)
    array = YAXArrays.Cube(filepath)
    geo_cube = GeoCube(array, latitude_name, longitude_name)
    return geo_cube
end

function GeoCube(data::Matrix, latitudes::AbstractVector, longitudes::AbstractVector)
    size(data)[1] == length(longitudes) || throw(ArgumentError("Matrix data must have the same number of rows than longitudes"))
    size(data)[2] == length(latitudes) || throw(ArgumentError("Matrix data must have the same number of columns than latitudes"))

    axlist = [
        RangeAxis("lon", longitudes),
        RangeAxis("lat", latitudes)
    ]
    geo_cube_arr = YAXArray(axlist, data)
    geo_cube = GeoCube(geo_cube_arr)
    return geo_cube
end

subsetcube(cube::GeoCube; kwargs...) = YAXArrays.Cubes.subsetcube(cube.data; kwargs...) |> GeoCube

function map_reduce_cells_to_geo(xout, xin::AbstractVector, cell_cube::CellCube, longitudes, latitudes)
    values_matrix = Matrix{eltype(xin)}(undef, length(longitudes), length(latitudes))

    for (lon_i, lon) in enumerate(longitudes)
        for (lat_i, lat) in enumerate(latitudes)
            cur_cell_id = get_cell_ids(cell_cube.grid, lat, lon)
            cur_cell_id in cell_cube.cell_ids || continue
            cur_cell_id_pos = findfirst(isequal(cur_cell_id), cell_cube.cell_ids)
            values_matrix[lon_i, lat_i] = cell_cube.data[cur_cell_id_pos]
        end
    end

    xout .= values_matrix
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

    # Expand spatial dimensions
    geo_array = mapCube(
        map_reduce_cells_to_geo,
        cell_cube.data,
        cell_cube,
        longitudes,
        latitudes,
        indims=InDims(:cell_id),
        outdims=OutDims(RangeAxis(:lon, longitudes), RangeAxis(:lat, latitudes))
    )
    cube = GeoCube(geo_array)
    return cube
end

function plot_map(geo_cube::GeoCube)
    axnames = [YAXArrays.Axes.axname(a) for a in geo_cube.data.axes]

    length(axnames) == 2 || throw(ArgumentError("Cube must have 2 axes lon and lat"))
    "lon" in axnames || throw(ArgumentError("Axis lon must be present"))
    "lat" in axnames || throw(ArgumentError("Axis lat must be present"))

    is_ordered = filter(x -> x in ["lon", "lat"], axnames) == ["lon", "lat"]
    data = is_ordered ? Matrix(geo_cube.data.data) : Matrix(geo_cube.data.data)'

    fig = Figure()
    ax = GeoAxis(fig[1, 1]; dest="+proj=wintri", coastlines=true)
    plt = surface!(ax, geo_cube.longitudes, geo_cube.latitudes, data; colormap=:viridis, shading=false)
    cb1 = Colorbar(fig[1, 2], plt; label="Value", height=Relative(0.5))
    return fig
end

function map_reduce_geo_to_cells(xout, xin; cell_ids_matrix, cell_ids::AbstractVector, aggregate_function::Function)
    cell_values = Vector{eltype(xin)}(undef, length(cell_ids))
    # allow for missing cell ids
    for (i, cell_id) in enumerate(cell_ids)
        cell_coords = findall(isequal(cell_id), cell_ids_matrix)
        if isempty(cell_coords)
            continue
        end
        cell_values[i] = xin[cell_coords] |> filter(!ismissing) |> aggregate_function
    end
    xout .= cell_values
end

"""
Import geographical data cube into a DGGS

Transforms a data cube with spatial index dimensions longitude and latitude
into a data cube with the cell id as a single spatial index dimension.
Re-gridding is done using the average value of all geographical coordinates belonging to a particular cell defined by the grid specification `grid_spec`.
"""
function CellCube(geo_cube::GeoCube, grid::AbstractGrid; aggregate_function::Function=mean)
    cell_ids_matrix = get_cell_ids(grid, geo_cube.latitudes, geo_cube.longitudes)
    cell_ids = cell_ids_matrix |> unique |> sort

    # Reduce spatial dimensions
    cell_array = mapCube(map_reduce_geo_to_cells, geo_cube.data,
        indims=InDims(:lat, :lon),
        outdims=OutDims(RangeAxis(:cell_id, cell_ids));
        cell_ids_matrix=cell_ids_matrix,
        cell_ids=cell_ids,
        aggregate_function=aggregate_function
    )
    cell_cube = CellCube(cell_array, grid)
    return cell_cube
end

function CellCube(data::AbstractVector, grid::AbstractGrid)
    length(data) == length(grid) || throw(ArgumentError("Vector data must have the same length as there are grid cells ($(length(grid)))"))

    axlist = [
        RangeAxis("cell_id", 1:length(grid)),
    ]
    cell_cube_arr = YAXArray(axlist, data)
    cell_cube = CellCube(cell_cube_arr, grid)
    return cell_cube
end

function plot_map(cell_cube::CellCube)
    axnames = [YAXArrays.Axes.axname(a) for a in cell_cube.data.axes]
    axnames == ["cell_id"] || throw(ArgumentError("Cube must have only cell_is as axis."))

    # raster cellcube, because geo coordinates needed for plotting
    geo_cube = GeoCube(cell_cube)
    plot_map(geo_cube)
end

function Base.getindex(cell_cube::CellCube; kwargs...)
    sub_cell_cube = YAXArrays.Cubes.subsetcube(cell_cube.data; kwargs...)
    CellCube(sub_cell_cube, cell_cube.grid)
end

Base.length(cell_cube::CellCube) = length(cell_cube.cell_ids)