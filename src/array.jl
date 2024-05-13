function DGGSArray(arr::YAXArray, name=:layer)
    haskey(arr.properties, "_DGGS") || error("Array is not in DGGS format")

    attrs = arr.properties
    level = attrs["_DGGS"]["level"]
    dggs = DGGSGridSystem(attrs["_DGGS"])

    DGGSArray(arr, attrs, name, level, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", arr::DGGSArray)
    println(io, "$(typeof(arr))")
    println(io, "DGGS: $(arr.dggs)")
    println(io, "Level: $(arr.level)")
    println(io, "Name: $(arr.name)")
end

function open_array(path::String)
    z = zopen(path)
    z isa ZArray || error("Path must point to a ZArray and not $(typeof(z))")
    data = zopen(path) |> YAXArray
    arr = YAXArray(data.axes, data, z.attrs)
    DGGSArray(arr)
end