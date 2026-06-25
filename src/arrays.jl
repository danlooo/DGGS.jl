function get_dggs_bbox(cells)
    cell = first(cells)
    resolution = cell.resolution

    # start with smallest possible bbox
    i_min = cell.i
    i_max = cell.i

    j_min, j_max = cell.j, cell.j
    n_min, n_max = cell.n, cell.n

    # extend bbox if needed
    for cell in cells
        if cell.i < i_min
            i_min = cell.i
        elseif cell.i > i_max
            i_max = cell.i
        end

        if cell.j < j_min
            j_min = cell.j
        elseif cell.j > j_max
            j_max = cell.j
        end

        if cell.n < n_min
            n_min = cell.n
        elseif cell.n > n_max
            n_max = cell.n
        end
    end

    return (
        Dim{:dggs_i}(i_min:i_max),
        Dim{:dggs_j}(j_min:j_max),
        Dim{:dggs_n}(n_min:n_max)
    )
end

"Infere max possible geo extent"
function get_geo_bbox(x::Union{DGGSArray,DGGSDataset})
    i_min, i_max = dims(x, :dggs_i).val.data |> x -> (first(x), last(x))
    j_min, j_max = dims(x, :dggs_j).val.data |> x -> (first(x), last(x))
    n_min, n_max = dims(x, :dggs_n).val.data |> x -> (first(x), last(x))

    dggs_corners = [
        Cell(i, j, n, x.resolution) for
        i in (i_min, i_max), j in (j_min, j_max), n in (n_min, n_max)
    ]
    geo_corners = to_geo.(dggs_corners)

    lon_min, lon_max = map(x -> x[1], geo_corners) |> x -> (minimum(x), maximum(x))
    lat_min, lat_max = map(x -> x[2], geo_corners) |> x -> (minimum(x), maximum(x))
    bbox = Extent(X=(lon_min, lon_max), Y=(lat_min, lat_max))
    return bbox
end

"Calculate actual geo extent"
function get_geo_bbox(geo_array::AbstractDimArray, crs::String; x_name=:X, y_name=:Y)
    # use default thread pool for lat/lon conversion
    wgs84_crs_geogcs = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]"

    x_min, x_max = dims(geo_array, x_name) |> extrema
    y_min, y_max = dims(geo_array, y_name) |> extrema

    if crs in [wgs84_crs_geogcs, "EPSG:4326"]
        ext = Extent(X=(x_min, x_max), Y=(y_min, y_max))
        return ext
    else
        trans = Proj.Transformation(crs, crs_geo)
        lat_min, lon_min = trans(x_min, y_min)
        lat_max, lon_max = trans(x_max, y_max)

        ext = Extent(X=(lon_min, lon_max), Y=(lat_min, lat_max))
        return ext
    end
end

function cells_to_coord_dict(cells::DimArray{Cell{Int64},2})
    cell_coords = Dict{Cell{Int64},Vector{CartesianIndex{2}}}()
    for cell_idx in CartesianIndices(cells)
        cell = cells[cell_idx]
        current_cells = get!(() -> CartesianIndex{2}[], cell_coords, cell)
        push!(current_cells, cell_idx)
    end
    return cell_coords
end


function to_dggs_array(
    geo_array::AbstractDimArray,
    cells,
    cell_coords,
    geo_bbox::Extent,
    agg_func::Function
    ;
    name=get_name(geo_array),
    out_eltype=Union{Missing,eltype(geo_array)},
    chunk_length=2^12,
    kwargs...
)
    resolution = first(cells).resolution
    dggsrs = "ISEA4D.Penta"

    # Create spatial dims
    spatial_dims = (Dim{:dggs_i}(0:(2*2^resolution-1)), Dim{:dggs_j}(0:(2^resolution-1)), Dim{:dggs_n}(0:4))

    # Create a TileArray directly instead of DGGSArray to avoid YAXArray wrapper overhead in setindex
    data = TileArray{out_eltype}(missing, length.(spatial_dims), (chunk_length, chunk_length, 1))

    # Pre-sized reusable buffer to avoid per-cell allocation when collecting non-missing values.
    # Most cells at high resolution map to ~1-4 pixels, so size 32 avoids reallocation in most cases.
    buf = Vector{eltype(geo_array)}(undef, 32)
    buf_len = 0

    # dims start at 0; +1 for 1-based Julia array indexing
    for (k, v) in cell_coords
        try
            # Collect non-missing values into pre-sized buffer (avoids geo_array[v] allocation)
            buf_len = 0
            for idx in v
                val = geo_array[idx]
                if val !== missing
                    buf_len += 1
                    if buf_len > length(buf)
                        resize!(buf, length(buf) * 2)
                    end
                    @inbounds buf[buf_len] = val
                end
            end
            buf_len == 0 && continue
            res = agg_func(@view buf[1:buf_len])
            data[k.i+1, k.j+1, k.n+1] = res
        catch
        end
    end

    return DGGSArray(
        data, spatial_dims, (), name, metadata(geo_array),
        resolution, dggsrs, geo_bbox
    )
end


function to_dggs_array(
    geo_array::AbstractDimArray, resolution::Integer, crs::String, agg_func::Function;
    x_name=:X, y_name=:Y, kwargs...
)
    x_dim = filter(x -> name(x) == x_name, dims(geo_array))
    y_dim = filter(x -> name(x) == y_name, dims(geo_array))
    isempty(x_dim) && error("X dimension (e.g. longitude) not found")
    isempty(y_dim) && error("Y dimension (e.g. latitude) not found")
    x_dim = only(x_dim)
    y_dim = only(y_dim)

    properties = metadata(geo_array)
    delete!(properties, "projection")

    cells = to_cell_array(x_dim, y_dim, resolution, crs)
    cell_coords = cells_to_coord_dict(cells)
    geo_bbox = get_geo_bbox(geo_array, crs)

    dggs_array = to_dggs_array(
        geo_array::AbstractDimArray,
        cells,
        cell_coords,
        geo_bbox::Extent,
        agg_func::Function
        ;
        name=get_name(geo_array),
        kwargs...
    )
    return dggs_array
end

function to_dggs_array(geo_array::AbstractDimArray, resolution::Integer, crs::String; x_name=:X, y_name=:Y, kwargs...)
    x_dim = filter(x -> name(x) == x_name, dims(geo_array))
    y_dim = filter(x -> name(x) == y_name, dims(geo_array))
    isempty(x_dim) && error("X dimension (e.g. longitude) not found")
    isempty(y_dim) && error("Y dimension (e.g. latitude) not found")
    x_dim = only(x_dim)
    y_dim = only(y_dim)

    properties = metadata(geo_array)

    cells = to_cell_array(x_dim, y_dim, resolution, crs)
    dggs_bbox = get_dggs_bbox(cells)
    geo_bbox = get_geo_bbox(geo_array, crs)

    dggs_array = to_dggs_array(geo_array, cells, dggs_bbox, geo_bbox; x_name=x_name, y_name=y_name, kwargs...)
    return dggs_array
end

function to_geo_array(dggs_array::DGGSArray, cells::AbstractDimArray; backend=:array, kwargs...)
    lon_dim = dims(cells, :X)
    lat_dim = dims(cells, :Y)

    # dggs_array may only contain parts of the world, having only parts of the dimension
    get_extent(i_dim) = dggs_array.dims[i_dim].val |> x -> (first(x), last(x))
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
                dggs_array[c][1]
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

function parse_bbox(bbox)
    if bbox isa Dict{String,Any}
        if haskey(bbox, "bounds")
            bbox = bbox["bounds"]
        end
        bbox = Extent(X=(bbox["X"]), Y=(bbox["Y"]))
    elseif bbox isa Extent
        # do nothing
    else
        bbox = Extent(bbox)
    end
    return bbox
end

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
    bbox = properties["dggs_bbox"] |> parse_bbox

    for k in ["dggs_resolution", "dggs_bbox", "dggs_dggsrs", "_FillValue", "fill_value", "missing_value"]
        delete!(properties, k)
    end

    arr_name = DD.name(array)
    if arr_name == DD.NoName()
        arr_name = get_name(array)
    end

    DGGSArray(
        array.data, dims(array), refdims(array), arr_name, properties,
        resolution, dggsrs, bbox
    )
end

function DGGSArray(resolution; chunk_length=2^12)
    spatial_dims = (Dim{:dggs_i}(0:(2*2^resolution-1)), Dim{:dggs_j}(0:(2^resolution-1)), Dim{:dggs_n}(0:4))
    data = TileArray{Union{Missing,Float64}}(missing, length.(spatial_dims), (chunk_length, chunk_length, 1))
    dggsrs = "ISEA4D.Penta"
    bbox = Extent(X=(-180, 180), Y=(-90, 90))
    a = DGGSArray(data, spatial_dims, (), DimensionalData.NoName(), DimensionalData.Dimensions.Lookups.NoMetadata(), resolution, dggsrs, bbox)
    return a
end

function YAXArrays.YAXArray(dggs_array::DGGSArray)
    properties = Dict{String,Any}(metadata(dggs_array))
    properties["dggs_resolution"] = dggs_array.resolution
    properties["dggs_dggsrs"] = dggs_array.dggsrs
    properties["dggs_bbox"] = NamedTuple(dggs_array.bbox)

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

Base.getindex(a::DGGSArray, c::Cell) = YAXArray(a)[dggs_i=At(c.i), dggs_j=At(c.j), dggs_n=At(c.n)]

# DGGSArrays are usually big. Like YAXArrays, avoid DiskArray to load everything in memory
Base.getindex(a::DGGSArray; i...) = view(a; i...)

Base.setindex!(a::DGGSArray, val, c::Cell) = YAXArray(a)[dggs_i=At(c.i), dggs_j=At(c.j), dggs_n=At(c.n)] = val


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

#
# Operations
#

"Determines if hyperrectangle a with dims i, j and n shares area with rectangle b"
function intersects(ai::UnitRange, bi::UnitRange, aj::UnitRange, bj::UnitRange, an::Integer, bn::Integer)
    if an != bn
        return false
    end
    overlap_i = first(ai) <= first(bi) <= last(ai) <= last(bi)
    overlap_j = first(aj) <= first(bj) <= last(aj) <= last(bj)
    overlap = overlap_i && overlap_j
    return overlap
end

"""
Crops the DGGSArray `a` to the dimensions of another DGGSArray `b`.
Resulting dimensions will be the intersection of those of `a` and `b`.
Returns a view into `a`.
"""
function crop(a::DGGSArray, b::DGGSArray)
    sel = Dict()

    spatial_dims = [:dggs_i, :dggs_j, :dggs_n]
    for dim in spatial_dims
        a_range = dims(a, dim) |> extrema |> x -> UnitRange(x[1], x[2])
        b_range = dims(b, dim) |> extrema |> x -> UnitRange(x[1], x[2])
        shared_range = intersect(a_range, b_range)
        sel[dim] = Between(shared_range.start, shared_range.stop)
    end

    # other dims, especially those that arn't UnitRange
    for dim in setdiff(name.(dims(a)), spatial_dims)
        shared = intersect(dims(a, dim), dims(b, dim))
        sel[dim] = At(shared)
    end
    return view(a; sel...)
end