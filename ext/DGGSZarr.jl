module DGGSZarr

using DGGS
using Zarr
using YAXArrays
using Infiltrator
using OrderedCollections
using Extents

# Currently, YAXArrays does not support saving the experimental nested DimTree type

function DGGS.save_dggs_pyramid(path::String, dggs_p::DGGSPyramid, args...; storetype=DirectoryStore)
    pyramid_attrs = Dict(
        "dggs_bbox" => dggs_p.bbox,
        "dggs_dggsrs" => dggs_p.dggsrs,
        "dggs_resolutions" => keys(dggs_p.data)
    )
    store = storetype(path, args...)
    group = zgroup(store; attrs=pyramid_attrs)
    for (resolution, dggs_ds) in dggs_p.data
        ds = Dataset(dggs_ds)
        savedataset(ds; path="$(path)/dggs_s$(resolution)")
    end
    return path
end

function DGGS.open_dggs_pyramid(path::String, args...; storetype=DirectoryStore)
    store = storetype(path, args...)
    group = zopen(store)
    dggsrs = get(group.attrs, "dggs_dggsrs", missing)
    bbox = get(group.attrs, "dggs_bbox", missing)
    resolutions = get(group.attrs, "dggs_resolutions", missing)

    ismissing(dggsrs) && error("DGGSRS not found in the pyramid metadata")
    ismissing(bbox) && error("Bounding box not found in the pyramid metadata")
    ismissing(resolutions) && error("Resolutions not found in the pyramid metadata")

    bbox = Extent(X=(bbox["X"][1], bbox["X"][2]), Y=(bbox["Y"][1], bbox["Y"][2]))

    data = OrderedDict{Int,DGGSDataset}()
    for resolution in resolutions
        ds_path = "$(path)/dggs_s$(resolution)"
        ds = open_dataset(ds_path)
        dggs_ds = DGGSDataset(ds)
        data[resolution] = dggs_ds
    end

    return DGGSPyramid(data, dggsrs, bbox)
end

end
