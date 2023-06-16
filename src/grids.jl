using NearestNeighbors
using YAXArrays
using Statistics
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

function Grid(grid_spec::GridSpec)
    data = grid_spec |> get_grid_data |> get_kd_tree
    grid = Grid(grid_spec, data)
    return grid
end

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
Get cell ids given geographic corrdinates
"""
function get_cell_ids(grid::Grid, lat_range::Union{AbstractVector,Number}, lon_range::Union{AbstractVector,Number})
    cell_ids = Matrix{Int}(undef, length(lat_range), length(lon_range))

    for (lat_i, lat) in enumerate(lat_range)
        for (lon_i, lon) in enumerate(lon_range)
            cell_ids[lat_i, lon_i] = NearestNeighbors.nn(grid.data, [lon, lat])[1]
        end
    end

    if size(cell_ids) == (1, 1)
        return cell_ids[1, 1]
    else
        return cell_ids
    end
end

"""
Convert cell id to geographic coordinate of cell center
"""
function get_geo_coords(grid::Grid, id::Int)
    tree_id = findfirst(x -> x == id, grid.data.indices)
    res = grid.data.data[tree_id]
    return (res[2], res[1])
end

function get_geo_coords(grid::Grid, id::AbstractVector)
    return get_geo_coords.(Ref(grid), id)
end

"""
Import geographical data cube into a DGGS

Transforms a data cube with spatial index dimensions longitude and latitude
into a data cube with the cell id as a single spatial index dimension.
Re-gridding is done using the average value of all geographical coordinates belonging to a particular cell defined by the grid specification `grid_spec`.
"""
function get_cell_cube(grid::Grid, geo_cube::YAXArray; latitude_name::String="lat", longitude_name::String="lon", aggregate_function::Function=mean)
    latitude_axis = getproperty(geo_cube, Symbol(latitude_name))
    longitude_axis = getproperty(geo_cube, Symbol(longitude_name))

    cell_ids = get_cell_ids(grid, latitude_axis, longitude_axis)
    cell_value_type = eltype(geo_cube)
    cell_values = Vector{Union{cell_value_type,Missing}}(missing, length(grid))

    for cell_id in unique(cell_ids)
        cell_coords = findall(isequal(cell_id), cell_ids)
        if isempty(cell_coords)
            continue
        end
        cell_values[cell_id] = aggregate_function(geo_cube.data'[cell_coords])
    end

    axlist = [RangeAxis("cell_id", range(1, length(grid)))]
    cell_cube = YAXArray(axlist, cell_values)
    return cell_cube
end

"""
Export cell data cube into a traditional geographical one

Transforms a data cube with one spatial index dimensions, i. e., the cell id,
into a traditional geographical data cube with two spatial index dimensions longitude and latitude.
Values are taken from the nearest cell.
"""
function get_geo_cube(grid::Grid, cell_cube::YAXArray)
    longitudes = -180:180
    latitudes = -90:90

    cell_value_type = eltype(cell_cube)
    regridded_matrix = Matrix{Union{cell_value_type,Missing}}(missing, length(longitudes), length(latitudes))

    for (lon_i, lon) in enumerate(longitudes)
        for (lat_i, lat) in enumerate(latitudes)
            cur_cell_id = get_cell_ids(grid, lat, lon)
            regridded_matrix[lon_i, lat_i] = cell_cube.data[cur_cell_id]
        end
    end

    axlist = [
        RangeAxis("lon", longitudes),
        RangeAxis("lat", latitudes)
    ]

    cube = YAXArray(axlist, regridded_matrix)
    return cube
end

"""
Create a grid system of different resolutions
"""
function create_grids(projection::String, aperture::Int, topology::String, n_resolutions::Int)
    grids = [Grid(projection, aperture, topology, resolution) for resolution in 0:n_resolutions-1]
    return grids
end

"""
Get id of parent cell
"""
function get_parent_cell_id(grids::Vector{Grid}, resolution::Int, cell_id::Int)
    if resolution == 1
        error("Lowest rosolutio can not have any further parents")
    end

    parent_resolution = resolution - 1
    geo_coord = get_geo_coords(grids[resolution], cell_id)
    parent_cell_id = get_cell_ids(grids[parent_resolution], geo_coord[1], geo_coord[2])
    return parent_cell_id
end

"""
Get ids of all child cells
"""
function get_children_cell_ids(grids::Vector{Grid}, resolution::Int, cell_id::Int)
    if resolution == length(grids)
        error("Highest resolution can not have any further children")
    end
    parent_cell_ids = get_parent_cell_id.(Ref(grids), resolution + 1, 1:length(grids[resolution+1].data.data))
    findall(x -> x == cell_id, parent_cell_ids)
end

"""
Get a cell data cube pyramid

Calculates a stack of cell data cubes with incrementally lower resolutions
based on the same data as provided by `cell_cube`.
Cell values are combined according to the provided `aggregate_function`.
"""
function get_cube_pyramid(grids::Vector{Grid}, cell_cube::YAXArray; aggregate_function::Function=mean)
    res = Vector{YAXArray}(undef, length(grids))
    res[length(grids)] = cell_cube

    # Calculate lower resolution based on the previous one
    for resolution in length(grids)-1:-1:1
        # parent: has higher resolution, used for combining
        # child: has lower resolution, to be calculated, stores the combined values
        parent_cell_cube = res[resolution+1]
        child_grid = grids[resolution]
        child_cell_vector = Vector{eltype(cell_cube)}(undef, length(child_grid))

        for cell_id in 1:length(child_grid)
            # downscaling by combining corresponding values from parent
            cell_ids = get_children_cell_ids(grids, resolution, cell_id)
            cell_values = [parent_cell_cube[cell_id=x].data for x in cell_ids]
            child_cell_vector[cell_id] = cell_values |> aggregate_function |> first
        end

        axlist = [
            RangeAxis("cell_id", range(1, length(child_grid)))
        ]
        res[resolution] = YAXArray(axlist, child_cell_vector)
    end
    return res
end

struct GridSystem
    grids::Vector{Grid}
    data::Vector{YAXArray}
    projection::String
    aperture::Int
    topology::String
    n_resolutions::Int
end

function GridSystem(geo_cube::YAXArray, projection::String, aperture::Int, topology::String, n_resolutions::Int; latitude_name::String="lat", longitude_name::String="lon")
    grids = create_grids(projection, aperture, topology, n_resolutions)
    finest_grid = grids[n_resolutions]
    finest_cell_cube = get_cell_cube(finest_grid, geo_cube; latitude_name=latitude_name, longitude_name=longitude_name)
    data = get_cube_pyramid(grids, finest_cell_cube)
    res = GridSystem(grids, data, projection, aperture, topology, n_resolutions)
    return res
end

# imported from https://github.com/JuliaDataCubes/YAXArrays.jl
cubesize(c::YAXArray{T}) where {T} = (sizeof(T)) * prod(map(length, caxes(c)))
cubesize(::YAXArray{T,0}) where {T} = sizeof(T)

function Base.show(io::IO, ::MIME"text/plain", grid_system::GridSystem)
    println(io, "Discrete Global Grid System")
    println(io, "Grid:\t$(grid_system.topology) topology, $(grid_system.projection) projection, aperture of $(grid_system.aperture)")
    println(io, "Cells:\t$(grid_system.n_resolutions) resolutions with up to $(grid_system.grids |> last |> length) cells")
    println(io, "Data:\tYAXArray of type $(typeof(grid_system.data[1].data)) with $(grid_system.data |> x -> map(cubesize, x) |> sum) bytes")
end


function get_geo_cube(dggs::GridSystem, resolution::Int=3)
    get_geo_cube(dggs.grids[resolution], dggs.data[resolution])
end

function get_cell_cube(dggs::GridSystem, resolution::Int=3)
    dggs.data[resolution]
end