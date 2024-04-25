struct ColorScale{T<:Real}
    schema::ColorScheme
    min_value::T
    max_value::T
end

struct Q2DI
    n
    i
    j
end

struct CellCube
    data::YAXArray
    level::Integer
end

struct GridSystem
    data::Dict{Int,CellCube} # some levels might be skipped

    function GridSystem(data)
        map(x -> x.data.axes .|> name, values(data)) |> allequal || error("Same dimensions must be used at all levels.")
        new(data)
    end
end