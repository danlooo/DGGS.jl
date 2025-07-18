module DGGSZarr

using DGGS
using Zarr
using YAXArrays
using FillArrays
using Infiltrator
using Extents
using OrderedCollections
using DimensionalData
import DimensionalData as DD
import Base.Threads: @threads

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

function DGGS.to_dggs_dataset(
    geo_ds::YAXArrays.Dataset,
    resolution::Integer,
    crs::String,
    store::Zarr.AbstractStore=DirectoryStore(tempname())
    ;
    outtype=Any,
    kwargs...
)
    # Create skeleton
    geo_bbox = DGGS.get_geo_bbox(geo_ds.cubes |> values |> first, crs
    )
    # Create first dataset skeletton
    ds_path = joinpath(store.folder, "dggs_s$(resolution)")
    ds_attrs = Dict(
        "dggs_bbox" => geo_bbox,
        "dggs_dggsrs" => "ISEA4D.Penta",
        "dggs_resolution" => resolution,
    )
    cells = DGGS.compute_cell_array(geo_ds.X, geo_ds.Y, resolution, crs)
    dggs_dims = DGGS.get_dggs_bbox(cells)

    skeleton_d = map(geo_ds.cubes |> keys |> collect) do k
        geo_arr = geo_ds.cubes[k]
        other_dims = filter(x -> !(name(x) in (:X, :Y)), dims(geo_arr))
        all_dims = isempty(other_dims) ? dggs_dims : vcat(dggs_dims, other_dims)
        (k => YAXArray(all_dims, Zeros(outtype, length.(all_dims)...), ds_attrs))
    end
    skeleton_ds = skeleton_d |> Dict |> x -> Dataset(; properties=ds_attrs, x...)
    savedataset(skeleton_ds; path=ds_path, driver=:zarr, skeleton=true, overwrite=true)

    @threads for array_key in collect(keys(geo_ds.cubes))
        geo_array = geo_ds.cubes[array_key]
        dggs_array = DGGS.to_dggs_array(
            geo_array, cells, dggs_dims, geo_bbox;
            outtype=outtype
        )
        z_arr = zopen(ds_path, "w")[String(array_key)]
        dggs_data = map(
            # Zarr does not support missing values directly
            x -> ismissing(x) ? z_arr.attrs["missing_value"] : x,
            dggs_array.data
        )
        # write array to disk
        z_arr[(Colon() for _ in 1:ndims(z_arr))...] = dggs_data
    end

    dggs_ds = open_dggs_dataset(ds_path; driver=:zarr)
    return dggs_ds
end

"Build pyramid, directly write to Zarr, optimized for parallel write of big datasets"
function DGGS.to_dggs_pyramid(
    geo_ds::YAXArrays.Dataset,
    resolution::Integer,
    crs::String,
    store::Zarr.AbstractStore
    ;
    pyramid_agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean,
    outtype=Any,
    outtype_sums=Any,
    kwargs...
)
    # Create pyramid skeleton
    geo_bbox = DGGS.get_geo_bbox(geo_ds.cubes |> values |> first, crs)
    pyramid_attrs = Dict(
        "dggs_bbox" => geo_bbox,
        "dggs_dggsrs" => "ISEA4D.Penta",
    )
    zgroup(store; attrs=pyramid_attrs)

    dggs_ds = DGGS.to_dggs_dataset(geo_ds, resolution, crs, store; outtype=outtype, outtype_sums=outtype_sums, kwargs...)
    # BUG: arrays were stored as dataset with subdir layer

    pyramid = DGGSDataset[]
    push!(pyramid, dggs_ds)
    current_dggs_ds = pyramid[end]
    for resolution in dggs_ds.resolution-1:-1:1
        coarser_arrays = []
        for array_key in keys(current_dggs_ds)
            dggs_array = getproperty(current_dggs_ds, array_key)
            coarser_dggs_array = DGGS.coarsen(
                dggs_array;
                path=joinpath(store.folder, "dggs_s$(resolution)"),
                pyramid_agg_func=pyramid_agg_func,
                kwargs...)
            push!(coarser_arrays, coarser_dggs_array)
        end
        coarser_ds = DGGSDataset(coarser_arrays...)
        current_dggs_ds = coarser_ds
        push!(pyramid, coarser_ds)
    end
    data = (pyramid |> reverse .|> x -> x.resolution => x) |> OrderedDict
    pyramid = DGGSPyramid(data, dggs_ds.dggsrs, dggs_ds.bbox)
    return pyramid
end
end
