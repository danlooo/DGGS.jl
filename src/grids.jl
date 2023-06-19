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
function get_geo_coords(grid::Grid, id::Int)
    tree_id = findfirst(x -> x == id, grid.data.indices)
    res = grid.data.data[tree_id]
    return (res[2], res[1])
end

"""
Convert cell id to geographic coordinate of cell center
"""
function get_geo_coords(grid::Grid, id::AbstractVector)
    return get_geo_coords.(Ref(grid), id)
end