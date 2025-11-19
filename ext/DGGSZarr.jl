module DGGSZarr

using DGGS
using Zarr
using YAXArrays
using Infiltrator
using Extents
using DimensionalData
import DimensionalData as DD
using DiskArrays
using FillArrays

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
    pyramid = DGGS.open_dggs_pyramid(group)
end

function DGGS.open_dggs_pyramid(group::ZGroup)
    dggsrs = get(group.attrs, "dggs_dggsrs", missing)
    bbox = get(group.attrs, "dggs_bbox", missing)

    ismissing(dggsrs) && error("DGGSRS not found in the pyramid metadata")
    ismissing(bbox) && error("Bounding box not found in the pyramid metadata")

    bbox = Extent(X=(bbox["X"][1], bbox["X"][2]), Y=(bbox["Y"][1], bbox["Y"][2]))

    dimtree = DimTree()
    groups = sort(group.groups, by=k -> group.groups[k].attrs["dggs_resolution"])
    for (k, v) in pairs(groups)
        ds = open_dataset(v; driver=:zarr) |> cache
        dggs_ds = DGGSDataset(ds)
        setproperty!(dimtree, Symbol(k), dggs_ds)
    end
    res = DGGSPyramid(dimtree, dggsrs, bbox)
    return res
end

function DGGS.init_global_dggs_dataset(geo_ds::Dataset, resolution, crs, path; x_dim_name=:X, y_dim_name=:Y, kwargs...)
    # extract spatial dimensions
    all_dims = []
    for (k, c) in geo_ds.cubes
        append!(all_dims, dims(c))
    end
    x_dim = filter(x -> name(x) == x_dim_name, all_dims)[1]
    y_dim = filter(x -> name(x) == y_dim_name, all_dims)[1]
    bbox = DGGS.get_geo_bbox(x_dim, y_dim, crs)

    properties = Dict(
        "dggs_resolution" => resolution,
        "dggs_dggsrs" => "ISEA4D.Penta",
        "dggs_bbox" => bbox
    )
    x_dim_name in keys(geo_ds.axes) || error("x_dim_name :$(x_dim_name) not found in geo_ds")
    y_dim_name in keys(geo_ds.axes) || error("y_dim_name :$(y_dim_name) not found in geo_ds")
    x_dim_name != y_dim_name || error("X and Y names must be different")

    arrays = Dict()
    for (key, geo_array) in pairs(geo_ds.cubes)
        if x_dim_name in name(geo_array.axes) && y_dim_name in name(geo_array.axes)
            spatial_dims = (Dim{:dggs_i}(0:2*2^resolution-1), Dim{:dggs_j}(0:2^resolution-1), Dim{:dggs_n}(0:4))
            non_spatial_dims = filter(x -> !(name(x) in [x_dim_name, y_dim_name]), geo_array.axes)
            dims = (spatial_dims..., non_spatial_dims...)

            chunk_size = Int(floor(2^resolution / 2))
            data = Zeros(eltype(geo_array), length.(dims))
            chunks = DiskArrays.GridChunks(data, (chunk_size, chunk_size, 1, fill(1, length(non_spatial_dims))...))
        else
            spatial_dims = ()
            non_spatial_dims = geo_array.axes
            dims = (spatial_dims..., non_spatial_dims...)
            data = Zeros(eltype(geo_array), length.(dims))
            chunks = DiskArrays.Unchunked()
        end

        yax_array = YAXArray(dims, data, properties, chunks=chunks)
        yax_array = rebuild(yax_array; name=key)
        arrays[key] = yax_array
    end

    ds = Dataset(; properties, arrays...)
    ds = savedataset(ds; path=path, skeleton=true, driver=:zarr, kwargs...)
    res = open_dataset(zopen(path, "w"); driver=:zarr) |> DGGSDataset
    return res
end
end