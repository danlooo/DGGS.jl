import YAXArrays: Cubes.formatbytes, Cubes.cubesize

struct Level{G<:AbstractGrid}
    data::CellCube
    grid::G
    level::Int
end

Base.length(level::Level) = length(level.data)
Base.getindex(level::Level, i...) = level.data[i...]
Base.eltype(level::Level) = eltype(level.data)

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
Get id of parent cell
"""
function get_parent_cell_id(grids::Vector{<:AbstractGrid}, level::Int, cell_id::Int)
    level > 1 || throw(ArgumentError("Lowest resolution can not have any further parents"))

    parent_level = level - 1
    geo_coord = get_geo_coords(grids[level], cell_id)
    parent_cell_id = get_cell_ids(grids[parent_level], geo_coord[1], geo_coord[2])
    return parent_cell_id
end

"""
Get ids of all child cells
"""
function get_children_cell_ids(grids::Vector{<:AbstractGrid}, level::Int, cell_id::Int)
    level < length(grids) || throw(ArgumentError("Highest level can not have any further children"))

    parent_cell_ids = get_parent_cell_id.(Ref(grids), level + 1, 1:length(grids[level+1].data.data))
    findall(x -> x == cell_id, parent_cell_ids)
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

    # Calculate lower level based on the previous one
    for level in length(grids)-1:-1:1
        # parent: has higher level, used for combining
        # child: has lower level, to be calculated, stores the combined values
        parent_cell_cube = res[level+1]
        current_grid = grids[level]
        child_cell_vector = Vector{eltype(cell_cube)}(undef, length(current_grid))

        for cell_id in 1:length(current_grid)
            # downscaling by combining corresponding values from parent
            cell_ids = get_children_cell_ids(grids, level, cell_id)
            values = parent_cell_cube[cell_ids]
            filtered_values = filter(x -> !ismissing(x), values)
            child_cell_vector[cell_id] = aggregate_function(filtered_values)
        end

        res[level] = CellCube(child_cell_vector, current_grid)
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