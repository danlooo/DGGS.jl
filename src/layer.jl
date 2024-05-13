
function DGGSLayer(data::YAXArrays.Dataset)
    haskey(data.properties, "_DGGS") || error("Dataset is not in DGGS format")

    layer = Dict{Symbol,DGGSArray}()
    for (k, c) in data.cubes
        arr = YAXArray(c.axes, c.data, union(c.properties, data.properties) |> Dict)
        layer[k] = DGGSArray(arr, k)
    end
    DGGSLayer(layer, data.properties)
end

function DGGSLayer(data::Dict{Symbol,DGGSArray}, attrs=Dict{String,Any}())
    bands = keys(data) |> collect
    level = data |> values |> first |> x -> x.level
    dggs = data |> values |> first |> x -> x.dggs
    DGGSLayer(data, attrs, bands, level, dggs)
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

Base.propertynames(l::DGGSLayer) = union(l.bands, (:data, :bands, :attrs))

function open_layer(path::String)
    z = zopen(path)
    z isa ZGroup || error("Path must point to a ZGoup and not $(typeof(z))")
    ds = open_dataset(z)
    DGGSLayer(ds)
end
