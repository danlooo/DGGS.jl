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
    println(io, "DGGS: $(arr.dggs)")
    println(io, "Level: $(arr.level)")
    println(io, "Name: $(arr.name)")
end

open_array(path::String) = error("Not implemented")