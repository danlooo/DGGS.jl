using NearestNeighbors
using DataFrames
using ArchGDAL


"""
A set of points defining a grid as a voronoi partition. 
Points must be geographical coordinates (lon, lat) and stored as a KDTree in the filed data.
"""
abstract type AbstractGrid end

struct Grid <: AbstractGrid
    data::KDTree

    function Grid(tree::KDTree)
        all([-180 <= x[1] <= 180 for x in tree.data]) ||
            throw(ArgumentError("Longitudes must be in range [-180, 180]"))
        all([-90 <= x[2] <= 90 for x in tree.data]) ||
            throw(ArgumentError("Latitudes must be in range [-90, 90]"))
        new(tree)
    end
end

"""
Create a Grid using coordinates of center points describing a voronoi partition.
center_points must have one point per column with 2 rows for longitude and latitude, respectiveley.
"""
function Grid(center_points::AbstractMatrix{<:Number})
    tree = KDTree(center_points)
    Grid(tree)
end

Base.length(grid::AbstractGrid) = Base.length(grid.data.data)

function Base.show(io::IO, ::MIME"text/plain", grid::AbstractGrid)
    print(io, "Grid with $(length(grid)) cells")
end

"""
Get cell ids of k nearest neighbors arround a cell
"""
function knn(grid::AbstractGrid, lat::Real, lon::Real, k::Integer)
    res = NearestNeighbors.knn(grid.data, [lat, lon], k)[1]
    return res
end

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
    projection::Union{Symbol,Nothing}
    aperture::Union{Int,Nothing}
    topology::Union{Symbol,Nothing}
    resolution::Int
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
function DgGrid(projection::Symbol, aperture::Union{Int,Nothing}, topology::Union{Symbol,Nothing}, resolution::Int)
    projection in Projections || throw(ArgumentError("projection :$projection must be one of $Projections"))
    aperture in Apertures || throw(ArgumentError("aperture $aperture must be one of $Apertures"))
    topology in Topologies || throw(ArgumentError("topology :$(topology) must be one of $Topologies"))

    grid_table = get_dggrid_grid_table(topology, projection, resolution)

    # cell center points encode grid tpopology (e.g. hexagon or square) implicitly
    # Fast average search in O(log n) and efficient in batch processing
    # KDTree defaults to Euklidean metric
    # However, should be faster than haversine and return same indices
    grid_tree = grid_table[:, [:lon, :lat]] |> Matrix |> transpose |> KDTree
    return DgGrid(grid_tree, :custom, projection, aperture, topology, resolution)
end

function DgGrid(preset::Symbol, resolution::Int)
    preset in keys(Presets) || throw(ArgumentError("Symbol $(preset) must be one of $(Presets)"))

    grid_table = get_dggrid_grid_table(preset, resolution)

    # cell center points encode grid tpopology (e.g. hexagon or square) implicitly
    # Fast average search in O(log n) and efficient in batch processing
    # KDTree defaults to Euklidean metric
    # However, should be faster than haversine and return same indices
    grid_tree = grid_table[:, [:lon, :lat]] |> Matrix |> transpose |> KDTree
    grid_props = Presets[preset]
    return DgGrid(grid_tree, preset, grid_props.projection, grid_props.aperture, grid_props.topology, resolution)
end

function Base.show(io::IO, ::MIME"text/plain", grid::DgGrid)
    print(io, "DgGrid $(grid.type) with $(grid.topology) topology, $(grid.projection) projection, aperture of $(grid.aperture), and $(length(grid)) cells at resolution $(grid.resolution)")
end

function get_cell_boundaries(grid::DgGrid)
    get_dggrid_cell_boundaries(grid.topology, grid.projection, grid.resolution)
end

create_toy_grid() = DgGrid(:isea, 4, :hexagon, 3)
