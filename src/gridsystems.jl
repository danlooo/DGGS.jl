struct Level
    cube::CellCube
    grid::AbstractGrid
end

struct GlobalGridSystem
    data::Vector{Level}
    projection::String
    aperture::Int
    topology::String
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

function GlobalGridSystem(geo_cube::GeoCube, n_levels::Int=5, projection=:isea::Symbol, aperture=4::Int, topology::Symbol=:hexagon)
    grids = [DgGrid(projection, aperture, topology, level) for level in 0:n_levels-1]
    finest_grid = grids[n_levels]
    finest_cell_cube = CellCube(geo_cube, finest_grid)
    cell_cubes = get_cube_pyramid(grids, finest_cell_cube)
end