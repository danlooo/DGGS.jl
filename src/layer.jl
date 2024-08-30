
function DGGSLayer(data::YAXArrays.Dataset)
    haskey(data.properties, "dggs_id") || error("Dataset is not in DGGS format")

    layer = Dict{Symbol,DGGSArray}()
    for (k, c) in data.cubes
        arr = YAXArray(c.axes, c.data, union(c.properties, data.properties) |> Dict)
        layer[k] = DGGSArray(arr, k)
    end
    DGGSLayer(layer, data.properties)
end

function DGGSLayer(data::Dict{Symbol,DGGSArray}, attrs=Dict{String,Any}())
    level = data |> values |> first |> x -> x.level
    dggs = data |> values |> first |> x -> x.dggs
    DGGSLayer(data, level, attrs, dggs)
end

function DGGSLayer(arr::DGGSArray)
    data = Dict{Symbol,DGGSArray}()
    data[:layer] = arr
    DGGSLayer(data, arr.level, arr.attrs, arr.dggs)
end

function DGGSLayer(arrs::Vector{DGGSArray{T,L}}) where {T,L}
    [x.id for x in arrs] |> allunique || error("IDs of arrays must be different")
    [x.level for x in arrs] |> allequal || error("Level of arrays must be the same")

    arr = first(arrs)
    data = Dict(x.id => x for x in arrs)
    DGGSLayer(data, arr.level, arr.attrs, arr.dggs)
end

DGGSLayer(::Vector{DGGSArray{T}}) where {T} = error("Level of arrays must be the same")

function Base.axes(l::DGGSLayer)
    axes = Vector()
    for arr in values(l.data)
        append!(axes, arr.data.axes)
    end
    unique!(axes)
    return axes
end

function show_arrays(io::IO, arrs::Vector{DGGSArray})
    printstyled(io, "Arrays:\n"; color=:white)

    for a in arrs
        print(io, "  ")
        Base.show(io, a)
        println(io, "")
    end
end

function Base.show(io::IO, ::MIME"text/plain", l::DGGSLayer)
    printstyled(io, typeof(l); color=:white)
    println(io, "")

    if "title" in keys(l.attrs)
        println(io, "Title:\t$(l.attrs["title"])")
    end
    println(io, "DGGS:\t$(l.dggs) at level $(l.level)")
    show_axes(io, axes(l))
    show_arrays(io, l.data |> values |> collect)
end

function Base.getproperty(l::DGGSLayer, v::Symbol)
    arrs = getfield(l, :data)
    if v in keys(arrs)
        return arrs[v]
    else
        return getfield(l, v)
    end
end

Base.getindex(l::DGGSLayer, id::Symbol) = l.data[id]

function Base.getindex(l::DGGSLayer, args...; kwargs...)
    id = get(kwargs, :id, nothing)
    isnothing(id) && error("Array id not provided.")

    center = get(kwargs, :center, nothing)
    lon = get(kwargs, :lon, nothing)
    lat = get(kwargs, :lat, nothing)
    radii = get(kwargs, :radii, nothing)

    args = filter(!isnothing, (center, lon, lat, radii))
    kwargs = filter(x -> !(x.first in [:id, :level, :center, :lat, :lon, :radii]), kwargs)

    if isempty(args) & isempty(kwargs)
        return l[id]
    else
        return getindex(l[id], args...; kwargs...)
    end

end
Base.propertynames(l::DGGSLayer) = union(l.data |> keys, (:data, :attrs))

"""
Transforms a `YAXArrays.Dataset` in geographic lat/lon ratser to a DGGSLayer at agiven layer
"""
function to_dggs_layer(
    geo_ds::Dataset,
    level::Integer;
    lon_name::Symbol=:lon,
    lat_name::Symbol=:lat,
    verbose::Bool=true,
    cell_ids::Union{AbstractMatrix,Nothing}=nothing,
    kwargs...)
    lon_dim = filter(x -> name(x) == lon_name, geo_ds.axes |> values |> collect)
    lat_dim = filter(x -> name(x) == lat_name, geo_ds.axes |> values |> collect)

    isempty(lon_dim) && error("Longitude dimension not found")
    isempty(lat_dim) && error("Latitude dimension not found")
    lon_dim = lon_dim[1]
    lat_dim = lat_dim[1]

    verbose && @info "Transform coordinates"
    cell_ids = isnothing(cell_ids) ? transform_points(lon_dim.val, lat_dim.val, level) : cell_ids

    data = Dict{Symbol,DGGSArray}()
    for (band, geo_arr) in geo_ds.cubes
        verbose && @info "Tranform band $band"
        data[band] = to_dggs_array(geo_arr, level; cell_ids=cell_ids, verbose=false, lon_name=lon_name, lat_name=lat_name, kwargs...)
    end
    properties = geo_ds.properties
    properties = union(properties, Q2DI_DGGS_PROPS) |> Dict
    properties["dggs_level"] = level

    DGGSLayer(data, level, properties, DGGSGridSystem(Q2DI_DGGS_PROPS))
end

to_dggs_layer(raster::AbstractDimArray, level::Integer, args...; kwargs...) = to_dggs_array(raster, level; kwargs...) |> DGGSLayer