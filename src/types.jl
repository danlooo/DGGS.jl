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

struct DGGSArray
    data::YAXArray
    level::Integer
end

struct DGGSArrayPyramid
    data::Dict{Int,DGGSArray} # some levels may be skipped

    function DGGSArrayPyramid(data)
        map(x -> x.data.axes .|> name, values(data)) |> allequal || error("Same dimensions must be used at all levels.")
        new(data)
    end
end

struct DGGSDataset
    data::YAXArrays.Dataset
end

struct DGGSDatasetPyramid
    data::Dict{Integer,DGGSDataset}
end

struct DGGSGridSystem
    name::String
    aperture::Int
    index::String
end

DGGSGridSystem(d::Dict{String,Any}) = DGGSGridSystem(d["name"], d["aperture"], d["index"])