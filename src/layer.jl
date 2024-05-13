
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
    dggs = data |> values |> first |> x -> x.dggs
    DGGSLayer(data, bands, level, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", l::DGGSLayer)
    println(io, "$(typeof(l))")
    println(io, "DGGS: $(l.dggs)")
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

open_layer(path::String) = error("Not implemented")
