Projections = ["ISEA", "FULLER"]
Topologies = ["HEXAGON", "TRIANGLE", "DIAMOND"]
GridPresets = ["SUPERFUND", "PLANETRISK", "ISEA4T", "ISEA4D", "ISEA3H", "ISEA4H", "ISEA7H", "ISEA43H", "FULLER4T", "FULLER4D", "FULLER3H", "FULLER4H", "FULLER7H", "FULLER43H"]
Apertures = [3, 4, 7]

struct Grid
    type::String
    projection::String
    aperture::Int
    topology::String
    resolution::Int
end

PresetParams = Dict(
    "FULLER7H" => Grid("FULLER7H", "FULLER", 7, "HEXAGON", 9),
    "ISEA4H" => Grid("ISEA4H", "ISEA", 4, "HEXAGON", 9)
)

function Grid(preset::String)
    if !(preset in GridPresets)
        throw(DomainError("preset must be an any of $(join(GridPresets, ","))"))
    end

    return PresetParams[preset]
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

    return Grid("CUSTOM", projection, aperture, topology, resolution)
end