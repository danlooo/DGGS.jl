using NearestNeighbors

Projections = ["ISEA", "FULLER"]
Topologies = ["HEXAGON", "TRIANGLE", "DIAMOND"]
GridPresets = ["SUPERFUND", "PLANETRISK", "ISEA4T", "ISEA4D", "ISEA3H", "ISEA4H", "ISEA7H", "ISEA43H", "FULLER4T", "FULLER4D", "FULLER3H", "FULLER4H", "FULLER7H", "FULLER43H"]
Apertures = [3, 4, 7]

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

PresetGridSpecs = Dict(
    "FULLER7H" => GridSpec("FULLER7H", "FULLER", 7, "HEXAGON", 9),
    "ISEA4H" => GridSpec("ISEA4H", "ISEA", 4, "HEXAGON", 9)
)

function Grid(preset::String)
    if !(preset in GridPresets)
        throw(DomainError("preset must be an any of $(join(GridPresets, ","))"))
    end

    spec = PresetGridSpecs[preset]
    data = generate_centers(spec)
    return Grid(spec, data)
end

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
    data = generate_centers(spec)
    return Grid(spec, data)
end

"Convert geographic corrdinates to cell id"
function cell_name(grid::Grid, lat::Real, lon::Real)
    if abs(lat) > 90
        throw(DomainError("Latitude argument lat must be within [-90, 90]"))
    end

    if abs(lon) > 180
        throw(DomainError("Longitude argument lon must be within [-180, 180]"))
    end

    NearestNeighbors.nn(grid.data, [lat, lon])[1]
end

"Convert cell id to geographic coordinate of cell center"
geo_coords(grid::Grid, ids::Vector{Int}) = grid.data.data[ids]