
function DGGSPyramid(data::Dict{Int,DGGSLayer})
    levels = data |> keys |> collect
    bands = data |> values |> first |> x -> x.bands
    dggs = data |> values |> first |> x -> x.dggs
    DGGSPyramid(data, levels, bands, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DGGSPyramid)
    println(io, "$(typeof(dggs))")
    println(io, "DGGS: $(dggs.dggs)")
    println(io, "Levels: $(dggs.levels)")
    println(io, "Bands: $(dggs.bands)")
end

Base.getindex(dggs::DGGSPyramid, i::Integer) = dggs.data[i]


function open_pyramid(path::String)
    root_group = zopen(path)
    haskey(root_group.attrs, "_DGGS") || error("Zarr store is not in DGGS format")

    pyramid = Dict{Int,DGGSLayer}()
    for level in root_group.attrs["_DGGS"]["levels"]
        layer_ds = open_dataset(root_group.groups["$level"])
        pyramid[level] = DGGSLayer(layer_ds)
    end

    return DGGSPyramid(pyramid)
end


function write_pyramid(path::String, dggs::DGGSPyramid)
    error("Not implemented")
end

function to_pyramid(data::AbstractDimArray)
    error("Not implemented")
end