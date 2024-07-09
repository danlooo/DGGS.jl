
function DGGSLayer(data::YAXArrays.Dataset)
    haskey(data.properties, "_DGGS") || error("Dataset is not in DGGS format")

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
    DGGSLayer(data, attrs, level, dggs)
end

function DGGSLayer(arr::DGGSArray)
    data = Dict{Symbol,DGGSArray}()
    data[:layer] = arr
    DGGSLayer(data, arr.attrs, arr.level, arr.dggs)
end

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

Base.getindex(l::DGGSLayer, array::Symbol) = l.data[array]
function Base.getindex(l::DGGSLayer; id::Symbol, kwargs...)
    if isempty(kwargs)
        return l[id]
    else
        return Base.getindex(l[id]; kwargs...)
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

    lon_dim = filter(x -> x isa X || name(x) == lon_name, collect(values(geo_ds.axes)))
    lat_dim = filter(x -> x isa Y || name(x) == lat_name, collect(values(geo_ds.axes)))

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
    DGGSLayer(data, geo_ds.properties, level, DGGSGridSystem(Q2DI_DGGS_PROPS))
end

to_dggs_layer(raster::AbstractDimArray, level::Integer; kw...) = to_dggs_array(raster, level, kw...) |> DGGSLayer