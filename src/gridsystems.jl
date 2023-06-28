import YAXArrays: Cubes.formatbytes, Cubes.cubesize

struct Level{G<:AbstractGrid}
    data::CellCube
    grid::G
    level::Int
end

Base.length(level::Level) = length(level.data)
Base.eltype(level::Level) = eltype(level.data)

function Base.getindex(level::Level; kwargs...)
    new_data = Base.getindex(level.data; kwargs...)
    result = Level(new_data, level.grid, level.level)
    return result
end

function Base.show(io::IO, ::MIME"text/plain", level::Level)
    println(io, "DGGS Level")
    println(io, "Cells:   level $(level.level) with $(length(level.data)) cells of type $(eltype(eltype(level.data)))")
    println(io, "Grid:    $(strip(repr("text/plain", level.grid)))")
    print(io, "Size:    $(formatbytes(cubesize(level.data)))")
end

plot_map(level::Level) = plot_map(level.data)

abstract type AbstractGlobalGridSystem <: AbstractVector{Level} end

struct GlobalGridSystem <: AbstractGlobalGridSystem
    data::Vector{Level}
    function GlobalGridSystem(levels)
        length(levels) >= 1 || throw(ArgumentError("Must provide at least one level"))
        new(levels)
    end
end

Base.length(dggs::AbstractGlobalGridSystem) = length(dggs.data)
Base.getindex(dggs::AbstractGlobalGridSystem, i...) = dggs.data[i...]
Base.lastindex(dggs::AbstractGlobalGridSystem) = lastindex(dggs.data)
Base.last(dggs::AbstractGlobalGridSystem) = last(dggs.data)
Base.size(dggs::AbstractGlobalGridSystem) = size(dggs.data)

function Base.show(io::IO, ::MIME"text/plain", dggs::AbstractGlobalGridSystem)
    println(io, "DGGS $(typeof(dggs))")
    println(io, "Cells: $(eltype(dggs.data))")
    println(io, "Grid:  $(strip(repr("text/plain", dggs[1].grid)))")
    print(io, "Levels:    $(length(dggs))")
end

struct DgGlobalGridSystem <: AbstractGlobalGridSystem
    data::Vector{Level}
    type::Symbol
    projection::Union{Symbol,Missing}
    aperture::Union{Int,Missing}
    topology::Union{Symbol,Missing}
end

"""
Get the cell ids of the children associated to each cell id of the parent level

Note that a child often has multiple parents.
The cells at different resolutions are only in rectangular pyramids perfectly nested.
Returns a Dict with parent cell id as keys and children cell ids as values
"""
function get_children_cell_ids(grids::Vector{<:AbstractGrid}, parent_level::Int, n_children::Int)
    res = Dict{Int,Vector{Int}}()
    parent_grid = grids[parent_level]
    child_grid = grids[parent_level+1]
    for parent_cell_id in 1:length(parent_grid)
        parent_cell_coord = get_geo_coords(parent_grid, parent_cell_id)
        child_cell_ids = DGGS.knn(child_grid, parent_cell_coord[1], parent_cell_coord[2], n_children)
        res[parent_cell_id] = child_cell_ids
    end
    return res
end

function reduce_cells_to_lower_resolution(xout, xin, child_cell_cube, child_cell_ids, parent_cell_ids, aggregate_function)
    coarsed_values = Vector{eltype(xin)}(undef, length(parent_cell_ids))

    for parent_cell_id in parent_cell_ids
        # allow non global cubes with missing cell ids
        selected_child_ids = child_cell_ids[parent_cell_id]
        selected_child_positions = findall(x -> x in selected_child_ids, child_cell_cube.cell_ids)

        coarsed_values[parent_cell_id] =
            xin[selected_child_positions] |>
            x -> filter(!ismissing, x) |>
                 aggregate_function
    end

    xout .= coarsed_values
end

"""
Get a cell data cube pyramid

Calculates a stack of cell data cubes with incrementally lower levels
based on the same data as provided by `cell_cube`.
Cell values are combined according to the provided `aggregate_function`.
"""
function get_cube_pyramid(grids::Vector{<:AbstractGrid}, cell_cube::CellCube; aggregate_function::Function=mean)
    res = Vector{CellCube}(undef, length(grids))
    res[length(grids)] = cell_cube

    # Calculate parent layers by aggregating previous child layer
    for parent_level in length(grids)-1:-1:1
        parent_grid = grids[parent_level]
        parent_cell_ids = 1:length(parent_grid)
        child_level = parent_level + 1
        child_cell_ids = get_children_cell_ids(grids, parent_level, 7)
        child_cell_cube = res[child_level]

        parent_cell_array = mapCube(
            reduce_cells_to_lower_resolution,
            child_cell_cube.data,
            child_cell_cube,
            child_cell_ids,
            parent_cell_ids,
            aggregate_function,
            indims=InDims(:cell_id),
            outdims=OutDims(RangeAxis(:cell_id, parent_cell_ids))
        )

        res[parent_level] = CellCube(parent_cell_array, parent_grid)
    end

    return res
end

function DgGlobalGridSystem(geo_cube::GeoCube, preset::Symbol, n_levels::Int=5)
    grids = [DgGrid(preset, level) for level in 0:n_levels-1]
    finest_grid = grids[n_levels]
    finest_cell_cube = CellCube(geo_cube, finest_grid)
    cell_cubes = get_cube_pyramid(grids, finest_cell_cube)
    levels = [Level(cell_cube, grid, level) for (cell_cube, grid, level) in zip(cell_cubes, grids, 1:length(grids))]
    dggs = DgGlobalGridSystem(levels, preset, missing, missing, missing)
    return dggs
end

function DgGlobalGridSystem(geo_cube::GeoCube, n_levels::Int=5, projection::Symbol=:isea, aperture::Int=4, topology::Symbol=:hexagon)
    grids = [DgGrid(projection, aperture, topology, level) for level in 0:n_levels-1]
    finest_grid = grids[n_levels]
    finest_cell_cube = CellCube(geo_cube, finest_grid)
    cell_cubes = get_cube_pyramid(grids, finest_cell_cube)
    levels = [Level(cell_cube, grid, level) for (cell_cube, grid, level) in zip(cell_cubes, grids, 1:length(grids))]
    dggs = DgGlobalGridSystem(levels, :custom, projection, aperture, topology)
    return dggs
end

function get_apertures(dggs::DgGlobalGridSystem)
    apertures = [level.grid.aperture for level in dggs]
    if length(unique(apertures)) == 1
        return apertures[1]
    else
        return apertures
    end
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DgGlobalGridSystem)
    println(io, "DGGS $(dggs.type) $(typeof(dggs))")
    println(io, "Cells:   $(length(dggs)) levels with up to $(dggs |> last |> length) cells of type $(dggs |> last |> eltype)")
    println(io, "Grid:    DgGrid with $(dggs.topology) topology, $(dggs.projection) projection, and aperture $(get_apertures(dggs))")
    print(io, "Size:    $(formatbytes(sum([cubesize(x.data) for x in dggs])))")
end