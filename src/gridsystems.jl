struct Level
    data::CellCube
    grid::AbstractGrid
end

Base.getindex(level::Level, i...) = level.data[i...]

function Base.show(io::IO, ::MIME"text/plain", level::Level)
    println(io, "DGGS Level")
    println(io, "Cells: $(length(level.data)) cells of type $(eltype(eltype(level.data)))")
    println(io, "Grid:  $(repr("text/plain", level.grid))")
end

function Base.show(io::IO, ::MIME"text/plain", dggs::AbstractGlobalGridSystem)
    println(io, "DGGS $(typeof(dggs))")
    println(io, "Element type: $(eltype(1))")
end

abstract type AbstractGlobalGridSystem end

struct GlobalGridSystem <: AbstractGlobalGridSystem
    data::Vector{Level}
end

Base.eltype(dggs::AbstractGlobalGridSystem) = eltype(dggs[1].data)
Base.length(dggs::AbstractGlobalGridSystem) = length(dggs.data)
Base.getindex(dggs::AbstractGlobalGridSystem, i...) = dggs.data[i...]
Base.lastindex(dggs::AbstractGlobalGridSystem) = last(dggs.data)

function Base.show(io::IO, ::MIME"text/plain", dggs::AbstractGlobalGridSystem)
    println(io, "DGGS $(typeof(dggs))")
    println(io, "Element type: $(eltype(dggs))")
    println(io, "Levels:       $(length(dggs))")
end

struct DgGlobalGridSystem <: AbstractGlobalGridSystem
    data::Vector{Level}
    projection::Symbol
    aperture::Int
    topology::Symbol
end

"""
Get id of parent cell
"""
function get_parent_cell_id(grids::Vector{<:AbstractGrid}, resolution::Int, cell_id::Int)
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
function get_children_cell_ids(grids::Vector{<:AbstractGrid}, resolution::Int, cell_id::Int)
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
function get_cube_pyramid(grids::Vector{<:AbstractGrid}, cell_cube::CellCube; aggregate_function::Function=mean)
    res = Vector{CellCube}(undef, length(grids))
    res[length(grids)] = cell_cube

    # Calculate lower resolution based on the previous one
    for resolution in length(grids)-1:-1:1
        # parent: has higher resolution, used for combining
        # child: has lower resolution, to be calculated, stores the combined values
        parent_cell_cube = res[resolution+1]
        current_grid = grids[resolution]
        child_cell_vector = Vector{eltype(cell_cube)}(undef, length(current_grid))

        for cell_id in 1:length(current_grid)
            # downscaling by combining corresponding values from parent
            cell_ids = get_children_cell_ids(grids, resolution, cell_id)
            child_cell_vector[cell_id] =
                parent_cell_cube[cell_ids] |>
                aggregate_function |>
                first
        end

        res[resolution] = CellCube(child_cell_vector, current_grid)
    end

    return res
end

function DgGlobalGridSystem(geo_cube::GeoCube, n_levels::Int=5, projection::Symbol=:isea, aperture::Int=4, topology::Symbol=:hexagon)
    grids = [DgGrid(projection, aperture, topology, level) for level in 0:n_levels-1]
    finest_grid = grids[n_levels]
    finest_cell_cube = CellCube(geo_cube, finest_grid)
    cell_cubes = get_cube_pyramid(grids, finest_cell_cube)
    levels = [Level(cell_cube, grid) for (cell_cube, grid) in zip(cell_cubes, grids)]
    println(levels)
    dggs = DgGlobalGridSystem(levels, projection, aperture, topology)
    return dggs
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DgGlobalGridSystem)
    println(io, "DGGS $(typeof(dggs))")
    println(io, "Element type: $(eltype(dggs))")
    println(io, "Levels:       $(length(dggs)) levels with up to $(length(last(dggs.data).data)) cells")
    println(io, "Grid:         DgGrid with $(dggs.topology) topology, $(dggs.projection) projection, and aperture of $(dggs.aperture)")
end