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