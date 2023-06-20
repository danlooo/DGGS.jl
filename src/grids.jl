using NearestNeighbors
using DataFrames
using ArchGDAL

abstract type AbstractGrid end

# Most basic grid. These fileds must be included in all subtypes
struct Grid <: AbstractGrid
    data::KDTree
    level::Int
end

Base.length(grid::AbstractGrid) = Base.length(grid.data.data)

"""
Get cell ids given geographic corrdinates
"""
function get_cell_ids(grid::AbstractGrid, lat_range::Union{AbstractVector,Number}, lon_range::Union{AbstractVector,Number})
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
function get_geo_coords(grid::AbstractGrid, id::Int)
    tree_id = findfirst(x -> x == id, grid.data.indices)
    res = grid.data.data[tree_id]
    return (res[2], res[1])
end

"""
Convert cell id to geographic coordinate of cell center
"""
function get_geo_coords(grid::AbstractGrid, id::AbstractVector)
    return get_geo_coords.(Ref(grid), id)
end

struct DgGrid <: AbstractGrid
    data::KDTree
    type::Symbol
    projection::Symbol
    aperture::Int
    topology::Symbol
    level::Int
end

function get_cell_centers(grid::AbstractGrid)
    # Using ArchGDAL directly results in segfaults and code would be more complex
    geometry = Vector{ArchGDAL.IGeometry}(undef, length(grid))
    for i in eachindex(grid.data.data)
        geometry[i] = ArchGDAL.createpoint(grid.data.data[i][1], grid.data.data[i][2])
    end
    df = DataFrame(geometry=geometry, cell_id=grid.data.indices)
    sort!(df, :cell_id)
    return df
end


"""
Create a grid using DGGRID parameters
"""
function DgGrid(projection::Symbol, aperture::Int, topology::Symbol, level::Int)
    projection in Projections ? true : error("projection :$projection must be one of $Projections")
    aperture in Apertures ? true : error("aperture $aperture must be one of $Apertures")
    topology in Topologies ? true : error("topology :$(topology) must be one of $Topologies")

    grid_table = get_dggrid_grid_table(topology, projection, level)

    # cell center points encode grid tpopology (e.g. hexagon or square) implicitly
    # Fast average search in O(log n) and efficient in batch processing
    # KDTree defaults to Euklidean metric
    # However, should be faster than haversine and return same indices
    grid_tree = grid_table[:, [:lon, :lat]] |> Matrix |> transpose |> KDTree
    return DgGrid(grid_tree, :custom, projection, aperture, topology, level)
end

function Base.show(io::IO, ::MIME"text/plain", grid::DgGrid)
    println(io, "DgGrid with $(grid.topology) topology, $(grid.projection) projection, aperture of $(grid.aperture), and $(length(grid)) cells")
end

function get_cell_boundaries(grid::DgGrid)
    get_dggrid_cell_boundaries(grid.topology, grid.projection, grid.level)
end

create_toy_grid() = DgGrid(:isea, 4, :hexagon, 3)
