function get_dggs_bbox(x_dim, y_dim, resolution, crs=DGGS.crs_geo)
    # get geo edges
    edge_points = vcat(
        map(y -> (first(x_dim), y), y_dim.val),
        map(y -> (last(x_dim), y), y_dim.val),
        map(x -> (x, first(y_dim)), x_dim.val),
        map(x -> (x, last(y_dim)), x_dim.val)
    )

    # cache trans for speed up
    trans = Proj.Transformation(crs, crs_isea; ctx=Proj.proj_context_create(), always_xy=true)
    cells = map(x -> to_cell(x[1], x[2], resolution, trans), edge_points)

    extents = []
    for n in 0:4
        cur_cells = filter(c -> c.n == n, cells)
        length(cur_cells) == 0 && continue
        i_min, i_max = map(x -> x.i, cur_cells) |> extrema
        j_min, j_max = map(x -> x.j, cur_cells) |> extrema
        extent = Dict(
            :dggs_i => i_min:i_max,
            :dggs_j => j_min:j_max,
            :dggs_n => n:n
        )
        push!(extents, extent)
    end
    return extents
end

function get_geo_bbox(x_dim, y_dim, crs)
    # get geo edges
    edge_points = vcat(
        map(y -> (first(x_dim), y), y_dim.val),
        map(y -> (last(x_dim), y), y_dim.val),
        map(x -> (x, first(y_dim)), x_dim.val),
        map(x -> (x, last(y_dim)), x_dim.val)
    )

    # cache trans for speed up
    trans = Proj.Transformation(crs, crs_geo; ctx=Proj.proj_context_create(), always_xy=true)
    geo_points = map(x -> trans(x[1], x[2]), edge_points)
    extent = Extent(
        X=map(x -> x[1], geo_points) |> extrema,
        Y=map(x -> x[2], geo_points) |> extrema
    )
    return extent
end

function get_chunks(resolution, x_dim, y_dim, crs; chunks=(dggs_i=4096, dggs_j=4096, dggs_n=1))
    # start with dummy global array
    spatial_dims = (Dim{:dggs_i}(0:2*2^resolution-1), Dim{:dggs_j}(0:2^resolution-1), Dim{:dggs_n}(0:4))
    data = Zeros(UInt8, length.(spatial_dims))
    yax_array = YAXArray(spatial_dims, data, Dict())
    yax_array = setchunks(yax_array, chunks)

    # calc extents and chunks
    dggs_bboxes = get_dggs_bbox(x_dim, y_dim, resolution, crs)

    used_chunks = filter(yax_array.chunks) do c
        c_i_min = c[1].start - 1
        c_i_max = c[1].stop - 1
        c_j_min = c[2].start - 1
        c_j_max = c[2].stop - 1
        c_n = c[3].start - 1

        overlap = true
        for r in dggs_bboxes
            if c_n != r[:dggs_n].start
                overlap = false
                continue
            end
            if c_i_max < r[:dggs_i].start || c_i_min > r[:dggs_i].stop
                overlap = false
                continue
            end
            if c_j_max < r[:dggs_j].start || c_j_min > r[:dggs_j].stop
                overlap = false
                continue
            end
        end
        overlap
    end
    return used_chunks
end

function to_dggs_array(
    geo_array::AbstractDimArray,
    resolution,
    crs,
    ;
    path=tempname() * ".dggs.zarr",
    name=get_name(geo_array),
    kwargs...
)
    geo_array = cache(geo_array)
    geo_ds = Dataset(layer=geo_array)
    dggs_array = DGGS.init_global_dggs_dataset(geo_ds, resolution, crs, path; kwargs...).layer

    x_dim = dims(geo_array, :X)
    y_dim = dims(geo_array, :Y)
    chunks = get_chunks(resolution, x_dim, y_dim, crs)

    # multi-threading without task migration making proj thread safe
    Threads.@threads :static for chunk in chunks
        # cache for other dims, e.g. time steps
        trans = Proj.Transformation(DGGS.crs_isea, crs, ctx=Proj.proj_context_create(), always_xy=true)
        geo_coords = map(Iterators.product(chunk...)) do (i, j, n)
            cell = Cell(i - 1, j - 1, n - 1, resolution)
            x, y = to_geo(cell, trans)
        end

        # TODO: loop over other dims

        data = map(x -> geo_array[X=Near(x[1]), Y=Near(x[2])][1], geo_coords)
        dggs_array[dggs_i=chunk[1], dggs_j=chunk[2], dggs_n=chunk[3]] = data
    end
    return dggs_array
end

function to_geo_array(dggs_array::DGGSArray, cells::AbstractDimArray; backend=:array, kwargs...)
    lon_dim = dims(cells, :X)
    lat_dim = dims(cells, :Y)

    # dggs_array may only contain parts of the world, having only parts of the dimension
    get_extent(i_dim) = dggs_array.dims[i_dim].val.data |> x -> (first(x), last(x))
    i_min, i_max = get_extent(1)
    j_min, j_max = get_extent(2)
    n_min, n_max = get_extent(3)

    geo_array = if backend == :array
        # in memory calculation
        # mapCube can write to disk but can not utilize the cache
        if dggs_array.data isa DiskArrayTools.CFDiskArray
            dggs_array = cache(dggs_array)
        end
        map(cells) do c
            try
                dggs_array[c]
            catch
                missing
            end
        end |> x -> YAXArray(dims(x), x.data, Dict())
    else
        mapCube(
            dggs_array,
            indims=InDims(
                :dggs_i,
                :dggs_j,
                :dggs_n
            ),
            outdims=OutDims(lon_dim, lat_dim, backend=backend),
            kwargs...
        ) do xout, xin
            for ci in CartesianIndices(xout)
                try
                    cell_ci = cells[ci] |> x -> CartesianIndex(x.i - i_min + 1, x.j - j_min + 1, x.n - n_min + 1)
                    xout[ci] = xin[cell_ci]
                catch
                    # not data available for this pixel
                end
            end
        end
    end

    return geo_array
end

function to_geo_array(dggs_array::DGGSArray, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    cells = to_cell_array(lon_dim, lat_dim, dggs_array.resolution)
    return to_geo_array(dggs_array::DGGSArray, cells; kwargs...)
end

function to_geo_array(dggs_array, lon_range::AbstractRange, lat_range::AbstractRange; kwargs...)
    lon_dim = X(lon_range)
    lat_dim = Y(lat_range)
    to_geo_array(dggs_array, lon_dim, lat_dim; kwargs...)
end

#
# DGGSArray features
#

function DGGSArray(array::AbstractDimArray, resolution::Integer, dggsrs::String="ISEA4D.Penta", bbox::Extent=Extent(X=(-180, 180), Y=(-90, 90)); name=DD.name(array), metadata=metadata(array))
    return DGGSArray(
        array.data, dims(array), refdims(array), name, metadata,
        resolution, dggsrs, bbox
    )
end

function DGGSArray(array::AbstractDimArray)
    properties = Dict{String,Any}(metadata(array))

    "dggs_resolution" in keys(properties) || error("Missing dggs_resolution in metadata")
    "dggs_dggsrs" in keys(properties) || error("Missing dggs_dggsrs in metadata")
    "dggs_bbox" in keys(properties) || error("Missing dggs_bbox in metadata")

    resolution = properties["dggs_resolution"] |> Int
    dggsrs = properties["dggs_dggsrs"] |> String
    bbox = properties["dggs_bbox"] |> x -> x isa Extent ? x : Extent(X=(x["X"][1], x["X"][2]), Y=(x["Y"][1], x["Y"][2]))

    delete!(properties, "dggs_resolution")
    delete!(properties, "dggs_dggsrs")
    delete!(properties, "dggs_bbox")

    arr_name = DD.name(array)
    if arr_name == DD.NoName()
        arr_name = get_name(array)
    end

    DGGSArray(
        array.data, dims(array), refdims(array), arr_name, properties,
        resolution, dggsrs, bbox
    )
end

function YAXArrays.YAXArray(dggs_array::DGGSArray)
    properties = Dict{String,Any}(metadata(dggs_array))
    properties["dggs_resolution"] = dggs_array.resolution
    properties["dggs_dggsrs"] = dggs_array.dggsrs
    properties["dggs_bbox"] = dggs_array.bbox

    return YAXArray(dims(dggs_array), dggs_array.data, properties)
end

"rebuild immutable objects with new field values. Part of any AbstractDimArray."
function DD.rebuild(
    dggs_array::DGGSArray, data::AbstractArray, dims::Tuple, refdims::Tuple, name, metadata
)
    DGGSArray(data, dims, refdims, name, metadata, dggs_array.resolution, dggs_array.dggsrs, dggs_array.bbox)
end

function get_name(array::AbstractDimArray)
    # as implemented in python xarray
    # uses CF conventions
    isempty(metadata(array)) && return DD.NoName()
    haskey(metadata(array), "long_name") && return metadata(array)["long_name"] |> Symbol
    haskey(metadata(array), "standard_name") && return metadata(array)["standard_name"] |> Symbol
    haskey(metadata(array), "name") && return metadata(array)["name"] |> Symbol
    return DD.NoName()
end

DD.label(dggs_array::DGGSArray) = string(DD.name(dggs_array))
DD.name(dggs_array::DGGSArray) = dggs_array.name

function non_spatial_dims(dggs_array::DGGSArray)
    spatial_dim_names = [:dggs_i, :dggs_j, :dggs_n]
    filter(x -> !(name(x) in spatial_dim_names), dggs_array.dims)
end

Base.getindex(a::DGGSArray, c::Cell) = a[dggs_i=At(c.i), dggs_j=At(c.j), dggs_n=At(c.n)]

#
# IO:: Serialization of DGGS Arrays
#

function open_dggs_array(file_path::String)
    ds = open_dataset(file_path) |> cache
    length(ds.cubes) == 1 || error("Path contains more than one Array")

    arr_name, arr = first(ds.cubes)
    return DGGSArray(arr)
end

function save_dggs_array(file_path::String, dggs_array::DGGSArray; kwargs...)
    ds = Dataset(; Dict(DD.name(dggs_array) => YAXArray(dggs_array))...)
    savedataset(ds; path=file_path, kwargs...)
end