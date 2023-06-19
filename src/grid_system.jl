struct Level
    cube::CellCube
    grid::Grid
end

struct GlobalGridSystem
    data::Vector{Level}
end

