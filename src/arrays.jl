function compute_cell_array(lon_dim, lat_dim, resolution)
    -180 <= minimum(lon_dim) <= maximum(lon_dim) <= 180 || error("Longitude must be within [-90,90]")
    -90 <= minimum(lat_dim) <= maximum(lat_dim) <= 90 || error("Latitude must be within [-90,90]")

    [(lon, lat) for lon in lon_dim, lat in lat_dim] |>
    x -> to_cell(x, resolution) |> x -> DimArray(; data=x, dims=(lon_dim, lat_dim))
end

function compute_cell_array(x_dim, y_dim, resolution, crs)
    # convert back to EPSG:4326, then do normal to_cell
    # x,y : coordinates in given crs
    # row,col: position in x and y dim vectors

    # use default thread pool for lat/lon conversion
    crs == "EPSG:4326" && return compute_cell_array(x_dim, y_dim, resolution)

    transformations = Channel{Proj.Transformation}(Inf)
    for _ in 1:Threads.nthreads()
        put!(transformations, Proj.Transformation(crs, crs_geo; ctx=Proj.proj_context_create()))
    end

    cells = Matrix(undef, length(x_dim), length(y_dim))

    Threads.@threads for i in CartesianIndices((1:length(x_dim), 1:length(y_dim)))
        row, col = i.I
        x, y = x_dim[row], y_dim[col]
        trans = take!(transformations)
        lat, lon = trans(x, y)
        cell = to_cell(lon, lat, resolution)
        cells[row, col] = cell
        put!(transformations, trans)
    end

    cell_array = DimArray(; data=cells, dims=(x_dim, y_dim))

    return cell_array
end

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
        cell.i < i_min && (i_min = cell.i)
        cell.i > i_max && (i_max = cell.i)

        cell.j < j_min && (j_min = cell.j)
        cell.j > j_max && (j_max = cell.j)

        cell.n < n_min && (n_min = cell.n)
        cell.n > n_max && (n_max = cell.n)
    end

    bbox = (
        Dim{:dggs_i}(i_min:i_max),
        Dim{:dggs_j}(j_min:j_max),
        Dim{:dggs_n}(n_min:n_max)
    )

    return bbox
end

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

function to_dggs_array(
    geo_array,
    cells;
    agg_func::Function=mean,
    outtype=Float64,
    backend=:array,
    path=tempname() * ".dggs.zarr",
    name=get_name(geo_array),
    kwargs...
)
    resolution = first(cells).resolution

    # get pixels to aggregate for each cell
    cell_coords = Dict{eltype(cells),Vector{CartesianIndex{2}}}()
    for cell_idx in CartesianIndices(cells)
        cell = cells[cell_idx]
        current_cells = get!(() -> CartesianIndex{2}[], cell_coords, cell)
        push!(current_cells, cell_idx)
    end

    dggs_bbox = get_dggs_bbox(keys(cell_coords))

    # re-grid
    res = mapCube(
        # mapCube can't find axes of other AbstractDimArrays e.g. Raster
        YAXArray(dims(geo_array), geo_array.data, metadata(geo_array));
        indims=InDims(dims(geo_array, :X), dims(geo_array, :Y)),
        outdims=OutDims(
            dggs_bbox...,
            outtype=outtype,
            backend=backend,
            path=path
        ), kwargs...) do xout, xin
        for ci in CartesianIndices(xout)
            i, j, n = ci.I
            try
                cell = Cell(dggs_bbox[1][i], dggs_bbox[2][j], dggs_bbox[3][n], resolution)
                cells = cell_coords[cell]
                res = agg_func(view(xin, cells))
                xout[i, j, n] = res
            catch
                # fill gap by averaging available neighbors
                xmin = clamp(i, 2, size(xout, 1) - 1) - 1
                ymin = clamp(j, 2, size(xout, 2) - 1) - 1
                res = filter(!ismissing, xout[xmin:xmin+2, ymin:ymin+2, n]) |> agg_func
                xout[i, j, n] = res
            end
        end
    end

    return DGGSArray(
        res.data, dims(res), refdims(res), name, metadata(geo_array),
        resolution, "ISEA4D.Penta"
    )
end

function to_dggs_array(geo_array, resolution, crs::AbstractString; x_name=:X, y_name=:Y, kwargs...)
    x_dim = filter(x -> name(x) == x_name, dims(geo_array))
    y_dim = filter(x -> name(x) == y_name, dims(geo_array))
    isempty(x_dim) && error("X dimension (e.g. longitude) not found")
    isempty(y_dim) && error("Y dimension (e.g. latitude) not found")
    x_dim = only(x_dim)
    y_dim = only(y_dim)

    properties = metadata(geo_array)
    delete!(properties, "projection")

    cells = compute_cell_array(x_dim, y_dim, resolution, crs)
    dggs_array = to_dggs_array(geo_array, cells)
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

    geo_array = mapCube(
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

    return geo_array
end

function to_geo_array(dggs_array::DGGSArray, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    cells = compute_cell_array(lon_dim, lat_dim, dggs_array.resolution)
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

function DGGSArray(array::AbstractDimArray, resolution::Integer, dggsrs::String="ISEA4D.Penta"; name=DD.name(array), metadata=metadata(array))
    return DGGSArray(
        array.data, dims(array), refdims(array), name, metadata,
        resolution, dggsrs
    )
end

function DGGSArray(array::AbstractDimArray)
    properties = Dict{String,Any}(metadata(array))

    "dggs_resolution" in keys(properties) || error("Missing dggs_resolution in metadata")
    "dggs_dggsrs" in keys(properties) || error("Missing dggs_dggsrs in metadata")

    resolution = properties["dggs_resolution"] |> Int
    dggsrs = properties["dggs_dggsrs"] |> String

    delete!(properties, "dggs_resolution")
    delete!(properties, "dggs_dggsrs")

    arr_name = DD.name(array)
    if arr_name == DD.NoName()
        arr_name = get_name(array)
    end

    DGGSArray(
        array.data, dims(array), refdims(array), arr_name, properties,
        resolution, dggsrs
    )
end

function YAXArrays.YAXArray(dggs_array::DGGSArray)
    properties = Dict{String,Any}(metadata(dggs_array))
    properties["dggs_resolution"] = dggs_array.resolution
    properties["dggs_dggsrs"] = dggs_array.dggsrs

    return YAXArray(dims(dggs_array), dggs_array.data, properties)
end

"rebuild immutable objects with new field values. Part of any AbstractDimArray."
function DD.rebuild(
    dggs_array::DGGSArray, data::AbstractArray, dims::Tuple, refdims::Tuple, name, metadata
)
    DGGSArray(data, dims, refdims, name, metadata, dggs_array.resolution, dggs_array.dggsrs)
end

function get_name(array::AbstractDimArray)
    # as implemented in python xarray
    # uses CF conventions
    isempty(array.properties) && return DD.NoName()
    haskey(array.properties, "long_name") && return array.properties["long_name"] |> Symbol
    haskey(array.properties, "standard_name") && return array.properties["standard_name"] |> Symbol
    haskey(array.properties, "name") && return array.properties["name"] |> Symbol
    return DD.NoName()
end

DD.label(dggs_array::DGGSArray) = string(DD.name(dggs_array))
DD.name(dggs_array::DGGSArray) = dggs_array.name

function non_spatial_dims(dggs_array::DGGSArray)
    spatial_dim_names = [:dggs_i, :dggs_j, :dggs_n]
    filter(x -> !(name(x) in spatial_dim_names), dggs_array.dims)
end

#
# IO:: Serialization of DGGS Arrays
#

function open_dggs_array(file_path::String)
    ds = open_dataset(file_path)
    length(ds.cubes) == 1 || error("Path contains more than one Array")

    arr_name, arr = first(ds.cubes)
    return DGGSArray(arr)
end

function save_dggs_array(file_path::String, dggs_array::DGGSArray; kwargs...)
    ds = Dataset(; Dict(DD.name(dggs_array) => YAXArray(dggs_array))...)
    savedataset(ds; path=file_path, kwargs...)
end