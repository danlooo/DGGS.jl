#
# DGGSArray
#

function DGGSArray(arr::YAXArray)
    haskey(arr.properties, "_DGGS") || error("Array is not in DGGS format")

    attrs = arr.properties
    name = attrs["name"]
    level = attrs["_DGGS"]["level"]
    dggs = DGGSGridSystem(attrs["_DGGS"])

    DGGSArray(arr, attrs, name, level, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", arr::DGGSArray)
    println(io, "$(typeof(arr))")
    println(io, "Level: $(arr.level)")
    println(io, "Name: $(arr.name)")
end

#
# DGGSLayer
#

function DGGSLayer(data::YAXArrays.Dataset)
    haskey(data.properties, "_DGGS") || error("Dataset is not in DGGS format")

    layer = Dict{Symbol,DGGSArray}()
    for (k, c) in data.cubes
        layer[k] = DGGSArray(c)
    end
    DGGSLayer(layer)
end

function DGGSLayer(data::Dict{Symbol,DGGSArray})
    bands = keys(data) |> collect
    level = data |> values |> first |> x -> x.level
    DGGSLayer(data, bands, level)
end

function Base.show(io::IO, ::MIME"text/plain", l::DGGSLayer)
    println(io, "$(typeof(l))")
    println(io, "Level: $(l.level)")
    println(io, "Bands: $(l.bands)")
end

function Base.getproperty(l::DGGSLayer, v::Symbol)
    if v in getfield(l, :bands) # prevent stack overflow
        return l.data[v]
    else
        return getfield(l, v)
    end
end

Base.propertynames(l::DGGSLayer) = union(l.bands, (:data, :bands))

#
# DGGSPyramid
#

function DGGSPyramid(data::Dict{Int,DGGSLayer})
    levels = data |> keys |> collect
    bands = data |> values |> first |> x -> x.bands
    DGGSPyramid(data, levels, bands)
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DGGSPyramid)
    println(io, "$(typeof(dggs))")
    println(io, "Levels: $(dggs.levels)")
    println(io, "Bands: $(dggs.bands)")
end

Base.getindex(dggs::DGGSPyramid, i::Integer) = dggs.data[i]