struct ColorScale{T<:Real}
    schema::ColorScheme
    min_value::T
    max_value::T
end

struct Q2DI{T<:Integer}
    n::UInt8
    i::T
    j::T
end

Q2DI(n::Integer, i::T, j::T) where {T<:Integer} = Q2DI(UInt8(n), T(i), T(j))

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