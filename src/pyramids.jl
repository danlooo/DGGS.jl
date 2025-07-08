function coarsen(dggs_array::DGGSArray)
    coarser_level = dggs_array.resolution - 1

    coarser_arr = mapCube(
        dggs_array;
        indims=InDims(:dggs_i, :dggs_j),
        outdims=OutDims(
            Dim{:dggs_i}(range(0; step=1, length=2 * 2^coarser_level)),
            Dim{:dggs_j}(range(0; step=1, length=2^coarser_level))
        )
    ) do xout, xin
        dggs_n = xin.indices[3]

        # TODO: implement
        xout .= 1
    end

    properties = Dict{String,Any}(metadata(dggs_array))
    properties["dggs_dggsrs"] = dggs_array.dggsrs
    properties["dggs_resolution"] = coarser_level
    properties["dggs_bbox"] = dggs_array.bbox

    coarser_dggs_arr = YAXArray(dims(coarser_arr), coarser_arr.data, properties) |> DGGSArray
    return coarser_dggs_arr
end

function coarsen(dggs_ds::DGGSDataset)
    coarser_arrays = []
    Threads.@threads for key in keys(dggs_ds)
        dggs_array = getproperty(dggs_ds, key)
        coarser_dggs_array = coarsen(dggs_array)
        push!(coarser_arrays, coarser_dggs_array)
    end
    res = DGGSDataset(coarser_arrays...)
    return res
end

function to_dggs_pyramid(dggs_ds::DGGSDataset)
    pyramid = DGGSDataset[]
    push!(pyramid, dggs_ds)
    for resolution in dggs_ds.resolution-1:-1:1
        current_dggs_ds = pyramid[end]
        coarser_ds = coarsen(current_dggs_ds)
        push!(pyramid, coarser_ds)
    end
    data = (pyramid |> reverse .|> x -> x.resolution => x) |> OrderedDict
    pyramid = DGGSPyramid(data, dggs_ds.dggsrs, dggs_ds.bbox)
    return pyramid
end

Base.show(io::IO, p::DGGSPyramid) = print(io, "DGGSPyramid $(p.dggsrs) with resolutions $(first(p.data).second.resolution):$(last(p.data).second.resolution)")

open_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")
save_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")