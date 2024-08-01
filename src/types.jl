struct Q2DI{T<:Integer}
    n::UInt8
    i::T
    j::T
end

Q2DI(n::Integer, i::T, j::T) where {T<:Integer} = Q2DI(UInt8(n), T(i), T(j))

defaultfillval(x::Q2DI) = nothing
Base.show(io::IO, i::Q2DI) = print(io, "Q2DI($(i.n),$(i.i),$(i.j))")

function Base.show(io::IO, ::MIME"text/plain", i::Q2DI{T}) where {T<:Integer}
    println(io, "Q2DI($(i.n), $(i.i), $(i.j))")
end

struct DGGSGridSystem
    id::String
    index::String
    polygon::String
    aperture::Integer
    projection::String
end

DGGSGridSystem(d::Dict{String,Any}) = DGGSGridSystem(d["id"], d["index"], d["polygon"], d["aperture"], d["projection"])

Base.show(io::IO, ::MIME"text/plain", dggs::DGGSGridSystem) = Base.show_default(io, dggs)

function Base.show(io::IO, dggs::DGGSGridSystem)
    polygons = Dict("hexagon" => "⬢")
    print(io, "$(dggs.id) $(get(polygons, dggs.polygon, "?"))")
end

struct DGGSArray{T,L}
    data::YAXArray
    attrs::Dict{String,Any}
    id::Symbol
    level::Integer
    dggs::DGGSGridSystem

    function DGGSArray(data, attrs, id, level, dggs)
        dims(data, :q2di_i) |> isnothing && error("Dimension :q2di_i not found")
        dims(data, :q2di_j) |> isnothing && error("Dimension :q2di_j not found")
        dims(data, :q2di_n) |> isnothing && error("Dimension :q2di_n not found")
        length(data.q2di_n) == 12 || error("Dimension :q2di_n must have length 12")
        log2(length(data.q2di_i)) % 1 == 0 || error("Dimension :q2di_i must have a length of a power of 2")
        log2(length(data.q2di_j)) % 1 == 0 || error("Dimension :q2di_j must have a length of a power of 2")
        length(data.q2di_i) == length(data.q2di_j) || error("Dimensions :q2di_i and :q2di_j must have the same length")

        new{eltype(data.data),level}(data, attrs, id, level, dggs)
    end
end

struct DGGSLayer{L}
    data::Dict{Symbol,DGGSArray}
    attrs::Dict{String,Any}
    level::Integer
    dggs::DGGSGridSystem

    function DGGSLayer(data, level, attrs=Dict{String,Any}(), dggs=DGGSGridSystem(Q2DI_DGGS_PROPS))
        level > 0 || error("Level must be positive")

        if !(map(x -> x.dggs.id, collect(values(data))) |> allequal)
            error("DGGS are different")
        end

        if !(map(x -> x.level, collect(values(data))) |> allequal)
            error("Levels are different")
        end

        if !(map(x -> x.dggs, collect(values(data))) |> allequal)
            error("Grid Systems are different")
        end

        new{level}(data, attrs, level, dggs)
    end
end

struct DGGSPyramid
    data::Dict{Int,DGGSLayer}
    attrs::Dict{String,Any}
    levels::Vector{Integer}
    dggs::DGGSGridSystem

    function DGGSPyramid(data, attrs, levels, dggs)
        if !(map(l -> l.data |> keys, collect(values(data))) |> allequal)
            error("Arrays are different")
        end

        if !(map(l -> l.dggs, collect(values(data))) |> allequal)
            error("Grid Systems are different")
        end

        new(data, attrs, levels, dggs)
    end
end


const Q2DI_DGGS_PROPS = Dict(
    "index" => "Q2DI",
    "aperture" => 4,
    "rotation_lon" => 11.25,
    "polyhedron" => "icosahedron",
    "id" => "DGGRID ISEA4H Q2DI",
    "radius" => 6371007.180918475,
    "polygon" => "hexagon",
    "rotation_lat" => 58.2825,
    "projection" => "isea",
    "rotation_azimuth" => 0
)