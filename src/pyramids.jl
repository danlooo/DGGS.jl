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