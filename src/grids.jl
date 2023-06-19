import NearestNeighbors: KDTree

abstract type AbstractGrid end

struct Grid <: AbstractGrid
    data::KDTree
end

Base.length(grid::AbstractGrid) = Base.length(grid.data.data)