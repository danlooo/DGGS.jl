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

function Base.show(io::IO, ::MIME"text/plain", i::Q2DI{T}) where {T<:Integer}
    println(io, "Q2DI($(i.n), $(i.i), $(i.j))")
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