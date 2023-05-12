@enum Projection ISEA FULLER
@enum Topology HEXAGON TRIANGLE DIAMOND
@enum GridPreset SUPERFUND PLANETRISK ISEA4T ISEA4D ISEA3H ISEA4H ISEA7H ISEA43H FULLER4T FULLER4D FULLER3H FULLER4H FULLER7H FULLER43H

struct Grid
    type::String
    projection::Projection
    aperture::Int
    topology::Topology
    resolution::Int
end

PresetParams = Dict(
    FULLER7H => Grid("FULLER7H", FULLER, 7, HEXAGON, 9),
    ISEA4H => Grid("ISEA4H", ISEA, 4, HEXAGON, 9)
)

function Grid(preset::GridPreset)
    return PresetParams[preset]
end

function Grid(projection::Projection, aperture::Int, topology::Topology, resolution::Int)
    return Grid("CUSTOM", projection, aperture, topology, resolution)
end