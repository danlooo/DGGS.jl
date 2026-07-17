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

"""
    _local_ranges(r, i, chunk_size)

Compute local indices within a chunk for writing data.
Handles boundary chunks that may be smaller than chunk_size.
"""
function _local_ranges(r, i, chunk_size)
    ntuple(length(r)) do d
        start_local = r[d].start - (i[d] - 1) * chunk_size[d]
        end_local = r[d].stop - (i[d] - 1) * chunk_size[d]
        start_local:end_local
    end
end

"""
    _write_tile_array_to_zarr!(disk_array, tile_array, chunks)

Write only present tiles from a TileArray to a Zarr array.
Missing tiles (chunks that are `missing`) are not written to disk.
"""
function _write_tile_array_to_zarr!(disk_array, tile_array, chunks)
    if !isnothing(chunks)
        disk_array = setchunks(disk_array, chunks)
    end

    for (r, i) in zip(DGGS.ranges(tile_array), findall(!ismissing, tile_array.data))
        chunk_data = tile_array.data[i]
        local_ranges = _local_ranges(r, i, tile_array.chunk_size)
        disk_array[r...] .= chunk_data[local_ranges...]
    end
end

"""
    save_dggs_array(file_path, dggs_array; chunks=nothing, kwargs...)

Save a DGGSArray to Zarr format, optimized for TileArrays.
Missing tiles are not written to disk, leveraging Zarr's sparse chunk storage.
Uses `fill_value=nothing` so chunks can be filled with any value.
"""
function DGGS.save_dggs_array(file_path::String, dggs_array::DGGSArray; chunks=nothing, kwargs...)
    yax_array = YAXArray(dggs_array)
    if !isnothing(chunks)
        yax_array = setchunks(yax_array, chunks)
    end

    # Check if underlying data is a TileArray
    if dggs_array.data isa TileArray
        # Create skeleton with fill_value=nothing
        ds = Dataset(; Dict(DD.name(dggs_array) => yax_array)...)
        disk_ds = savedataset(ds; path=file_path, skeleton=true, driver=:zarr, fill_value=nothing, kwargs...)

        # Write only present tiles
        _write_tile_array_to_zarr!(disk_ds[DD.name(dggs_array)], dggs_array.data, chunks)
    else
        # Fallback for non-TileArray data
        ds = Dataset(; Dict(DD.name(dggs_array) => yax_array)...)
        savedataset(ds; path=file_path, driver=:zarr, fill_value=nothing, kwargs...)
    end
    return file_path
end

"""
    save_dggs_dataset(file_path, dggs_ds; chunks=(dggs_i=4096, dggs_j=4096, dggs_n=1), kwargs...)

Save a DGGSDataset to Zarr format, optimized for TileArrays.
Missing tiles are not written to disk, leveraging Zarr's sparse chunk storage.
Uses `fill_value=nothing` so chunks can be filled with any value.
"""
function DGGS.save_dggs_dataset(file_path::String, dggs_ds::DGGSDataset; chunks=(dggs_i=4096, dggs_j=4096, dggs_n=1), kwargs...)
    yax_ds = Dataset(dggs_ds)
    if !isnothing(chunks)
        yax_ds = setchunks(yax_ds, chunks)
    end

    # Check if any underlying data is a TileArray
    has_tile_array = any(k -> getproperty(dggs_ds, k).data isa TileArray, keys(dggs_ds))


    if has_tile_array
        # Create skeleton with fill_value=nothing
        disk_ds = savedataset(yax_ds; path=file_path, skeleton=true, driver=:zarr, fill_value=nothing, kwargs...)

        # Write only present tiles for each variable
        for key in keys(dggs_ds)
            if getproperty(dggs_ds, key).data isa TileArray
                _write_tile_array_to_zarr!(disk_ds[key], getproperty(dggs_ds, key).data, chunks)
            end
        end
    else
        # Fallback for non-TileArray data
        savedataset(yax_ds; path=file_path, driver=:zarr, fill_value=nothing, kwargs...)
    end
    return file_path
end

"""
    save_dggs_pyramid(path, dggs_p; storetype=DirectoryStore, chunks=(dggs_i=4096, dggs_j=4096, dggs_n=1), kwargs...)

Save a DGGSPyramid to Zarr format, optimized for TileArrays.
Missing tiles are not written to disk, leveraging Zarr's sparse chunk storage.
Uses `fill_value=nothing` so chunks can be filled with any value.
"""
function DGGS.save_dggs_pyramid(path::String, dggs_p::DGGSPyramid, args...; storetype=DirectoryStore, chunks=(dggs_i=4096, dggs_j=4096, dggs_n=1), kwargs...)
    pyramid_attrs = Dict(
        "dggs_bbox" => dggs_p.bbox,
        "dggs_dggsrs" => dggs_p.dggsrs,
    )
    store = storetype(path, args...)
    group = zgroup(store; attrs=pyramid_attrs)

    for key in keys(dggs_p.branches)
        dggs_ds = getproperty(dggs_p, key) |> x -> x isa DGGSArray ? DGGSDataset(x) : x
        DGGS.save_dggs_dataset("$(path)/$(key)", dggs_ds; chunks=chunks, kwargs...)
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

    bbox = DGGS.parse_bbox(bbox)

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


"""
Need to init globally:
- allows parallel read and write 
- ij chunks often don't overlap eith neoighboring chunks on different dggs_n quad
- empty chunks are not stored on disk
"""
function DGGS.init_global_dggs_dataset(
    geo_ds::Dataset, resolution, path;
    bbox=(X=(-180, 180), Y=(-90, 90)),
    x_dim_name=:X, y_dim_name=:Y, chunks=(dggs_i=4096, dggs_j=4096, dggs_n=1), kwargs...
)
    # extract spatial dimensions
    all_dims = []
    for (k, c) in geo_ds.cubes
        append!(all_dims, dims(c))
    end

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
        is_spatial = x_dim_name in name(geo_array.axes) && y_dim_name in name(geo_array.axes)
        if is_spatial
            spatial_dims = (Dim{:dggs_i}(0:(2*2^resolution-1)), Dim{:dggs_j}(0:(2^resolution-1)), Dim{:dggs_n}(0:4))
            non_spatial_dims = filter(x -> !(name(x) in [x_dim_name, y_dim_name]), geo_array.axes)
            dims = (spatial_dims..., non_spatial_dims...)
        else
            spatial_dims = ()
            non_spatial_dims = geo_array.axes
            dims = (spatial_dims..., non_spatial_dims...)
        end

        data = Zeros(Union{Missing,eltype(geo_array)}, length.(dims))
        yax_array = YAXArray(dims, data, properties)
        yax_array = rebuild(yax_array; name=key)
        yax_array = setchunks(yax_array, chunks)

        arrays[key] = yax_array
    end

    ds = Dataset(; properties, arrays...)
    ds = savedataset(ds; path=path, skeleton=true, driver=:zarr, kwargs...)
    res = open_dataset(zopen(path, "w"); driver=:zarr) |> DGGSDataset
    return res
end

end
