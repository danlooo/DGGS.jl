function DGGSArray(arr::YAXArray, name=:layer)
    haskey(arr.properties, "_DGGS") || error("Array is not in DGGS format")

    attrs = arr.properties
    level = attrs["_DGGS"]["level"]
    dggs = DGGSGridSystem(attrs["_DGGS"])

    DGGSArray(arr, attrs, name, level, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", arr::DGGSArray)
    println(io, "$(typeof(arr))")
    println(io, "DGGS: $(arr.dggs)")
    println(io, "Level: $(arr.level)")
    println(io, "Name: $(arr.name)")
end

function open_array(path::String)
    z = zopen(path)
    z isa ZArray || error("Path must point to a ZArray and not $(typeof(z))")
    data = zopen(path) |> YAXArray
    arr = YAXArray(data.axes, data, z.attrs)
    DGGSArray(arr)
end

"Apply function f after filtering of missing and NAN values"
function filter_null(f)
    x -> x |> filter(!ismissing) |> filter(!isnan) |> f
end

function map_geo_to_cell_array(xout, xin, cell_ids_indexlist, agg_func)
    for (cell_id, cell_indices) in cell_ids_indexlist
        # xout is not a YAXArray anymore
        xout[cell_id.i+1, cell_id.j+1, cell_id.n+1] = agg_func(view(xin, cell_indices))
    end
end

function to_array(
    raster::AbstractDimArray,
    level::Integer;
    lon_name::Symbol=:lon,
    lat_name::Symbol=:lat,
    cell_ids::Union{AbstractMatrix,Nothing}=nothing,
    agg_func::Function=filter_null(mean),
    verbose::Bool=true
)
    level > 0 || error("Level must be positive")

    lon_dim = filter(x -> x isa X || name(x) == lon_name, dims(raster))
    lat_dim = filter(x -> x isa Y || name(x) == lat_name, dims(raster))

    isempty(lon_dim) && error("Longitude dimension not found")
    isempty(lat_dim) && error("Latitude dimension not found")
    lon_dim = lon_dim[1]
    lat_dim = lat_dim[1]

    lon_dim.val.order == DimensionalData.Dimensions.LookupArrays.ForwardOrdered() || error("Longitude must be sorted in forward ascending oder")
    lat_dim.val.order == DimensionalData.Dimensions.LookupArrays.ForwardOrdered() || error("Latitude must be sorted in forward ascending oder")
    -180 <= minimum(lon_dim) <= maximum(lon_dim) <= 180 || error("$(name(lon_dim)) must be within [-180, 180]")
    -90 <= minimum(lat_dim) <= maximum(lat_dim) <= 90 || error("$(name(lon_dim)) must be within [-90, 90]")

    verbose && @info "Step 1/2: Transform coordinates"
    cell_ids_mat = isnothing(cell_ids) ? transform_points(lon_dim.val, lat_dim.val, level) : cell_ids

    cell_ids_indexlist = Dict()
    for c in 1:size(cell_ids_mat, 2)
        for r in 1:size(cell_ids_mat, 1)
            cell_id = cell_ids_mat[r, c]
            if cell_id in keys(cell_ids_indexlist)
                # multiple pixels per cell
                push!(cell_ids_indexlist[cell_id], CartesianIndex(r, c))
            else
                cell_ids_indexlist[cell_id] = [CartesianIndex(r, c)]
            end
        end
    end

    verbose && @info "Step 2/2: Re-grid the data"
    cell_array = mapCube(
        map_geo_to_cell_array,
        # mapCube can't find axes of other AbstractDimArrays e.g. Raster
        YAXArray(dims(raster), raster.data),
        cell_ids_indexlist,
        agg_func,
        indims=InDims(lon_dim, lat_dim),
        outdims=OutDims(
            Dim{:q2di_i}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_j}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_n}(0:11)
        ),
        showprog=true
    )
    props = raster |> metadata |> typeof == Dict{String,Any} ? deepcopy(metadata(raster)) : Dict{String,Any}()
    props["_DGGS"] = deepcopy(Q2DI_DGGS_PROPS)
    props["_DGGS"]["level"] = level
    cell_array = YAXArray(cell_array.axes, cell_array.data, props)
    DGGSArray(cell_array)
end