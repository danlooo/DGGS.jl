using NearestNeighbors
using YAXArrays
using Interpolations
using Statistics
using LinearAlgebra
using DataFrames

"""
Specification of a grid as defined by DGGRID
"""
struct GridSpec
    type::String
    projection::String
    aperture::Int
    topology::String
    resolution::Int
end

"""
A grid as defined by DGGRID with the actual grid points
"""
struct Grid
    spec::GridSpec
    data::KDTree
end

function Base.show(io::IO, ::MIME"text/plain", grid::Grid)
    println(io, "DGGS Grid with $(grid.spec.topology) topology, $(grid.spec.projection) projection, apterture of $(grid.spec.aperture), and $(length(grid.data.data)) cells")
end

"""
DGGRID grid presets
"""
PresetGridSpecs = Dict(
    "FULLER7H" => GridSpec("FULLER7H", "FULLER", 7, "HEXAGON", 9),
    "ISEA4H" => GridSpec("ISEA4H", "ISEA", 4, "HEXAGON", 9)
)

"""
Create a grid using DGGRID preset
"""
function Grid(preset::String)
    if !(preset in GridPresets)
        throw(DomainError("preset must be an any of $(join(GridPresets, ","))"))
    end

    spec = PresetGridSpecs[preset]
    data = spec |> get_grid_data |> get_kd_tree
    return Grid(spec, data)
end

Base.length(grid::Grid) = Base.length(grid.data.data)

Grid() = Grid("ISEA4H")

"""
Create a grid using DGGRID parameters
"""
function Grid(projection::String, aperture::Int, topology::String, resolution::Int)
    if !(projection in Projections)
        throw(DomainError("Argument projection must be an any of $(join(Projections, ","))"))
    end

    if !(aperture in Apertures)
        throw(DomainError("Argument aperture must be an any of $(join(Apertures, ","))"))
    end

    if !(topology in Topologies)
        throw(DomainError("Argument topology must be an any of $(join(Topologies, ","))"))
    end

    spec = GridSpec("CUSTOM", projection, aperture, topology, resolution)
    data = spec |> get_grid_data |> get_kd_tree
    return Grid(spec, data)
end

"""
represent cell table as kd-tree of center points
"""
function get_kd_tree(df::DataFrame; longitude_col=:lon, latitude_col=:lat)
    # cell center points encode grid tpopology (e.g. hexagon or square) implicitly
    # Fast average search in O(log n) and efficient in batch processing
    # KDTree defaults to Euklidean metric
    # However, should be faster than haversine and return same indices
    kd_tree = df[:, [longitude_col, latitude_col]] |> Matrix |> transpose |> KDTree
    return kd_tree
end

create_toy_grid() = Grid("ISEA", 4, "HEXAGON", 3)

"""
Convert geographic corrdinates to cell id
"""
function get_cell_ids(grid::Grid, lat::Real, lon::Real)
    if abs(lat) > 90
        throw(DomainError("Latitude argument lat must be within [-90, 90]"))
    end

    if abs(lon) > 180
        throw(DomainError("Longitude argument lon must be within [-180, 180]"))
    end

    NearestNeighbors.nn(grid.data, [lon, lat])[1]
end

"""
Convert cell id to geographic coordinate of cell center
"""
function get_geo_coords(grid::Grid, id::Int)
    tree_id = findfirst(x -> x == id, grid.data.indices)
    res = grid.data.data[tree_id]
    return (res[2], res[1])
end

"""
Import geographical data cube into a DGGS

Transforms a data cube with spatial index dimensions longitude and latitude
into a data cube with the cell id as a single spatial index dimension.
Re-gridding is done using the average value of all geographical coordinates belonging to a particular cell defined by the grid specification `grid_spec`.
"""
function get_cell_cube(grid_spec::GridSpec, geo_cube::YAXArray, latitude_name="lat", longitude_name="lon")
    latitude_axis = getproperty(geo_cube, Symbol(latitude_name))
    longitude_axis = getproperty(geo_cube, Symbol(longitude_name))
    cell_ids = DGGS.get_cell_ids(grid_spec, latitude_axis, longitude_axis)

    # binary matrix mapping geographic coordinates to cell ids
    geo_cell_mapping_matrix = cell_ids' .== unique(cell_ids)
    geo_cube_vector = geo_cube[:, :, 1] |> vec
    replace!(geo_cube_vector, missing => 0) # ignore missing values in average

    # calculate new values using a weighted average of all points 
    # Account for different number of points per cell in the division
    cell_cube_matrix = geo_cell_mapping_matrix * geo_cube_vector ./ sum(geo_cell_mapping_matrix, dims=2)
    cell_cube_vector = cell_cube_matrix[:, 1]

    axlist = [CategoricalAxis("cell_id", unique(cell_ids))]
    cell_cube = YAXArray(axlist, cell_cube_vector)
    return cell_cube
end

get_cell_cube(grid::Grid, geo_cube::YAXArray, latitude_name, longitude_name) = get_cell_cube(grid.spec, geo_cube::YAXArray, latitude_name, longitude_name)

"""
Export cell data cube into a traditional geographical one

Transforms a data cube with one spatial index dimensions, i. e., the cell id,
into a traditional geographical data cube with two spatial index dimensions longitude and latitude.
Re-gridding is done by mapping the cell value into the center point of its cell.
Interpolation of values is performed, because grid cells on the same meridian may have different latitudes.
"""
function get_geo_cube(grid_spec::GridSpec, cell_cube::YAXArray)
    df = get_grid_data(grid_spec)
    df.value = cell_cube.data
    sort!(df, [:lon, :lat])
    itp = linear_interpolation((sort(df.lon), sort(df.lat)), Diagonal(df.value))

    itr_lon = minimum(df.lon):2:maximum(df.lon)
    itr_lat = minimum(df.lat):2:maximum(df.lat)

    regridded_matrix = Array{Float64}(undef, length(itr_lon), length(itr_lat))
    for (lon_i, lon_val) in enumerate(itr_lon)
        for (lat_i, lat_val) in enumerate(itr_lat)
            regridded_matrix[lon_i, lat_i] = itp(lon_val, lat_val)
        end
    end

    axlist = [
        RangeAxis("lon", itr_lon),
        RangeAxis("lat", itr_lat)
    ]

    cube = YAXArray(axlist, regridded_matrix)
    return cube
end

get_geo_cube(grid::Grid, cell_cube::YAXArray) = get_geo_cube(grid.spec, cell_cube)