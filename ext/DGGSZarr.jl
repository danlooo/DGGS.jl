module DGGSZarr

using DGGS
using Zarr
using YAXArrays
using Infiltrator
using Extents
using DimensionalData
import DimensionalData as DD

# Currently, YAXArrays does not support saving the experimental nested DimTree type

function DGGS.save_dggs_pyramid(path::String, dggs_p::DGGSPyramid, args...; storetype=DirectoryStore)
    pyramid_attrs = Dict(
        "dggs_bbox" => dggs_p.bbox,
        "dggs_dggsrs" => dggs_p.dggsrs,
    )
    store = storetype(path, args...)
    group = zgroup(store; attrs=pyramid_attrs)
    for key in keys(dggs_p.branches)
        dggs_ds = getproperty(dggs_p, key) |> x -> x isa DGGSArray ? DGGSDataset(x) : x
        ds = Dataset(dggs_ds)
        savedataset(ds; path="$(path)/$(key)", driver=:zarr)
    end
    return path
end

function DGGS.open_dggs_pyramid(path::String, args...; storetype=DirectoryStore)
    store = storetype(path, args...)
    group = zopen(store)

    dggsrs = get(group.attrs, "dggs_dggsrs", missing)
    bbox = get(group.attrs, "dggs_bbox", missing)

    ismissing(dggsrs) && error("DGGSRS not found in the pyramid metadata")
    ismissing(bbox) && error("Bounding box not found in the pyramid metadata")

    bbox = Extent(X=(bbox["X"][1], bbox["X"][2]), Y=(bbox["Y"][1], bbox["Y"][2]))

    dimtree = DimTree()
    groups = sort(group.groups, by=k -> group.groups[k].attrs["dggs_resolution"])
    for (k, v) in pairs(groups)
        ds = open_dataset(v; driver=:zarr)
        dggs_ds = DGGSDataset(ds)
        setproperty!(dimtree, Symbol(k), dggs_ds)
    end
    res = DGGSPyramid(dimtree, dggsrs, bbox)
    return res
end

end
