struct Q2DI{T<:Integer}
    n::UInt8
    i::T
    j::T
end

Q2DI(n::Integer, i::T, j::T) where {T<:Integer} = Q2DI(UInt8(n), T(i), T(j))

function Base.show(io::IO, ::MIME"text/plain", i::Q2DI{T}) where {T<:Integer}
    println(io, "Q2DI($(i.n), $(i.i), $(i.j))")
end

struct DGGSGridSystem
    name::String
    index::String
end

DGGSGridSystem(d::Dict{String,Any}) = DGGSGridSystem(d["name"], d["index"])

struct DGGSArray
    data::YAXArray
    attrs::Dict{String,Any}
    name::String
    level::Integer
    dggs::DGGSGridSystem
end

struct DGGSLayer
    data::Dict{Symbol,DGGSArray}
    bands::Vector{Symbol}
    level::Integer
    dggs::DGGSGridSystem

    function DGGSLayer(data, bands, level, dggs)
        if !(map(x -> x.dggs.name, collect(values(data))) |> allequal)
            error("DGGS are different")
        end

        if !(map(x -> x.level, collect(values(data))) |> allequal)
            error("Levels are different")
        end

        if !(map(x -> x.dggs, collect(values(data))) |> allequal)
            error("Grid Systems are different")
        end

        new(data, bands, level, dggs)
    end
end

struct DGGSPyramid
    data::Dict{Int,DGGSLayer}
    levels::Vector{Integer}
    bands::Vector{Symbol}
    dggs::DGGSGridSystem

    function DGGSPyramid(data, levels, bands, dggs)
        if !(map(l -> l.bands, collect(values(data))) |> allequal)
            error("Bands are different")
        end

        if !(map(l -> l.dggs, collect(values(data))) |> allequal)
            error("Grid Systems are different")
        end

        new(data, levels, bands, dggs)
    end
end


const Q2DI_DGGS_PROPS = Dict(
    "index" => "Q2DI",
    "aperture" => 4,
    "rotation_lon" => 11.25,
    "polyhedron" => "icosahedron",
    "name" => "DGGRID ISEA4H Q2DI",
    "radius" => 6371007.180918475,
    "polygon" => "hexagon",
    "rotation_lat" => 58.2825,
    "projection" => "+isea",
    "rotation_azimuth" => 0
)