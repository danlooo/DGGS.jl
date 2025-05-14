function compute_cell_array(lon_dim, lat_dim, resolution)
    -180 <= minimum(lon_dim) <= maximum(lon_dim) <= 180 || error("Longitude must be within [-90,90]")
    -90 <= minimum(lat_dim) <= maximum(lat_dim) <= 90 || error("Latitude must be within [-90,90]")

    [(lon, lat) for lon in lon_dim, lat in lat_dim] |>
    x -> to_cell(x, resolution)
end

function compute_cell_array(x_dim, y_dim, resolution, crs)
    # convert back to EPSG:4326, then do normal to_cell
    # x,y : coordinates in given crs
    # row,col: position in x and y dim vectors

    transformations = Channel{Proj.Transformation}(Inf)
    for _ in 1:Threads.nthreads()
        put!(transformations, Proj.Transformation(crs, crs_geo; ctx=Proj.proj_context_create()))
    end

    cells = Matrix(undef, length(x_dim), length(y_dim))

    Threads.@threads for i in CartesianIndices((1:length(x_dim), 1:length(y_dim)))
        row, col = i.I
        x, y = x_dim[row], y_dim[col]
        trans = take!(transformations)
        lat, lon = trans(y, x)
        cell = to_cell(lon, lat, resolution)
        cells[row, col] = cell
        put!(transformations, trans)
    end

    return cells
end

function to_dggs_array(geo_array, resolution, crs::AbstractString="EPSG:4326"; agg_func::Function=mean, outtype=Float64, path=tempname() * ".dggs.zarr", lon_name=:X, lat_name=:Y, kwargs...)
    x_dim = filter(x -> name(x) == lon_name, dims(geo_array))
    y_dim = filter(x -> name(x) == lat_name, dims(geo_array))
    isempty(x_dim) && error("Longitude dimension not found")
    isempty(y_dim) && error("Latitude dimension not found")
    x_dim = only(x_dim)
    y_dim = only(y_dim)

    cells = if crs == "EPSG:4326"
        compute_cell_array(x_dim, y_dim, resolution)
    else
        compute_cell_array(x_dim, y_dim, resolution, crs)
    end

    # get pixels to aggregate for each cell
    cell_coords = Dict{eltype(cells),Vector{CartesianIndex{2}}}()
    for cell_idx in CartesianIndices(cells)
        cell = cells[cell_idx]
        current_cells = get!(() -> CartesianIndex{2}[], cell_coords, cell)
        push!(current_cells, cell_idx)
    end

    # re-grid
    res = mapCube(
        # mapCube can't find axes of other AbstractDimArrays e.g. Raster
        YAXArray(dims(geo_array), geo_array.data, metadata(geo_array));
        indims=InDims(x_dim, y_dim),
        outdims=OutDims(
            Dim{:dggs_i}(0:(2*2^resolution-1)),
            Dim{:dggs_j}(0:(2^resolution-1)),
            Dim{:dggs_n}(0:4),
            outtype=outtype,
            path=path
        ), kwargs...) do xout, xin
        for ci in CartesianIndices(xout)
            i, j, n = ci.I
            try
                cell = Cell(i - 1, j - 1, n - 1, resolution)
                res = agg_func(view(xin, cell_coords[cell]))
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
        res.data, dims(res), refdims(res), get_name(geo_array), metadata(geo_array),
        resolution, "ISEA4D.Penta"
    )
end

function to_geo_array(dggs_array::DGGSArray, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    cells = compute_cell_array(lon_dim, lat_dim, dggs_array.resolution)
    geo_array = mapCube(
        dggs_array,
        indims=InDims(
            :dggs_i,
            :dggs_j,
            :dggs_n
        ),
        outdims=OutDims(lon_dim, lat_dim),
        kwargs...
    ) do xout, xin
        xout .= map(x -> xin[x.i+1, x.j+1, x.n+1], cells)
    end

    return geo_array
end

function to_geo_array(dggs_array, lon_range::AbstractRange, lat_range::AbstractRange; kwargs...)
    lon_dim = X(lon_range)
    lat_dim = Y(lat_range)
    to_geo_array(dggs_array, lon_dim, lat_dim; kwargs...)
end

#
# DGGSArray features
#

function DGGSArray(array::AbstractDimArray, resolution::Integer, dggsrs::String)
    return DGGSArray(
        array.data, dims(array), refdims(array), name(array), metadata(array),
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


function Base.show(io::IO, mime::MIME"text/plain", dggs_array::DGGSArray)
    println(io, "DGGSArray{", eltype(dggs_array), "}")
    println(io, "DGGS: ", dggs_array.dggsrs, " at resolution ", dggs_array.resolution,
        " (", @sprintf("%.1e", prod(size(dggs_array.data))), " cells)")
    println(io, "Based on: ", join(size(dggs_array.data), "x"), " ", typeof(dggs_array.data).name.name)

    if length(dggs_array.dims) > 3
        println(io, "Additional dimensions:")
        for dim in non_spatial_dims(dggs_array)
            print(io, "   ")
            DD.Dimensions.print_dimname(io, dim)
            print(io, " $(minimum(dim):step(dim):maximum(dim))")
        end
        println(io, "")
    else
        println(io, "Additional dimensions: none")
    end

    if length(dggs_array.metadata) > 0
        println(io, "Meta data:")
        for (key, value) in dggs_array.metadata
            println(io, "   $key: $value")
        end
    else
        println(io, "No meta data")
    end
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