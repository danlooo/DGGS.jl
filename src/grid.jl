using NearestNeighbors

struct GridSpec
    type::String
    projection::String
    aperture::Int
    topology::String
    resolution::Int
end

struct Grid
    spec::GridSpec
    data::Any
end

function Base.show(io::IO, ::MIME"text/plain", grid::Grid)
    println(io, "DGGS Grid with $(grid.spec.topology) topology, $(grid.spec.projection) projection, apterture of $(grid.spec.aperture), and $(length(grid.data.data)) cells")
end

PresetGridSpecs = Dict(
    "FULLER7H" => GridSpec("FULLER7H", "FULLER", 7, "HEXAGON", 9),
    "ISEA4H" => GridSpec("ISEA4H", "ISEA", 4, "HEXAGON", 9)
)

function Grid(preset::String)
    if !(preset in GridPresets)
        throw(DomainError("preset must be an any of $(join(GridPresets, ","))"))
    end

    spec = PresetGridSpecs[preset]
    data = get_grid_data(spec)
    return Grid(spec, data)
end

Base.length(grid::Grid) = Base.length(grid.data.data)

Grid() = Grid("ISEA4H")

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
    data = get_grid_data(spec)
    return Grid(spec, data)
end

create_toy_grid() = Grid("ISEA", 4, "HEXAGON", 3)

"Convert geographic corrdinates to cell id"
function get_cell_name(grid::Grid, lat::Real, lon::Real)
    if abs(lat) > 90
        throw(DomainError("Latitude argument lat must be within [-90, 90]"))
    end

    if abs(lon) > 180
        throw(DomainError("Longitude argument lon must be within [-180, 180]"))
    end

    NearestNeighbors.nn(grid.data, [lat, lon])[1]
end

"Convert cell id to geographic coordinate of cell center"
get_geo_coords(grid::Grid, ids::Vector{Int}) = grid.data.data[ids]
get_geo_coords(grid::Grid, id::Int) = grid.data.data[id]