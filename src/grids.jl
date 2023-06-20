import NearestNeighbors: KDTree

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


"""
Create a grid using DGGRID parameters
"""
function DgGrid(projection::Symbol, aperture::Int, topology::Symbol, resolution::Int)
    projection in Projections ? true : error("projection :$projection must be one of $Projections")
    aperture in Apertures ? true : error("aperture $aperture must be one of $Apertures")
    topology in Topologies ? true : error("topology :$(topology) must be one of $Topologies")

    grid_table = get_dggrid_grid_table(:custom, topology, projection, resolution)

    # cell center points encode grid tpopology (e.g. hexagon or square) implicitly
    # Fast average search in O(log n) and efficient in batch processing
    # KDTree defaults to Euklidean metric
    # However, should be faster than haversine and return same indices
    grid_tree = grid_table[:, [:lon, :lat]] |> Matrix |> transpose |> KDTree
    return DgGrid(grid_tree, :custom, projection, aperture, topology, resolution)
end

function Base.show(io::IO, ::MIME"text/plain", grid::DgGrid)
    println(io, "DgGrid with $(grid.topology) topology, $(grid.projection) projection, aperture of $(grid.aperture), and $(length(grid)) cells")
end