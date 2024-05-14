
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
    bands = keys(data) |> collect
    level = data |> values |> first |> x -> x.level
    dggs = data |> values |> first |> x -> x.dggs
    DGGSLayer(data, attrs, bands, level, dggs)
end

function DGGSLayer(arr::DGGSArray)
    bands = [:layer]
    data = Dict{Symbol,DGGSArray}()
    data[:layer] = arr
    DGGSLayer(data, arr.attrs, bands, arr.level, arr.dggs)
end

function Base.show(io::IO, ::MIME"text/plain", l::DGGSLayer)
    println(io, "$(typeof(l))")
    println(io, "DGGS: $(l.dggs)")
    println(io, "Level: $(l.level)")
    println(io, "Bands: $(l.bands)")
end

Base.getindex(l::DGGSLayer, band::Symbol) = l.data[band]

function Base.getproperty(l::DGGSLayer, v::Symbol)
    if v in getfield(l, :bands) # prevent stack overflow
        return l.data[v]
    else
        return getfield(l, v)
    end
end

Base.propertynames(l::DGGSLayer) = union(l.bands, (:data, :bands, :attrs))

function open_dggs_layer(path::String)
    z = zopen(path)
    z isa ZGroup || error("Path must point to a ZGoup and not $(typeof(z))")
    ds = open_dataset(z)
    DGGSLayer(ds)
end

"""
Transforms a `YAXArrays.Dataset` in geographic lat/lon ratser to a DGGSLayer at agiven layer
"""
function to_dggs_layer(geo_ds::Dataset, level::Integer; lon_name=:lon, lat_name=:lat)
    level > 0 || error("Level must be positive")

    lon_dim = filter(x -> x isa X || name(x) == lon_name, collect(values(geo_ds.axes)))
    lat_dim = filter(x -> x isa Y || name(x) == lat_name, collect(values(geo_ds.axes)))

    isempty(lon_dim) && error("Longitude dimension not found")
    isempty(lat_dim) && error("Latitude dimension not found")
    lon_dim = lon_dim[1]
    lat_dim = lat_dim[1]

    cell_ids = transform_points(lon_dim.val, lat_dim.val, level)
    data = Dict{Symbol,DGGSArray}()
    for (band, geo_arr) in geo_ds.cubes
        data[band] = to_dggs_array(geo_arr, level; cell_ids=cell_ids)
    end
    bands = geo_ds.cubes |> keys |> collect
    DGGSLayer(data, geo_ds.properties, bands, level, DGGSGridSystem(Q2DI_DGGS_PROPS))
end

to_dggs_layer(raster::AbstractDimArray, level::Integer; kw...) = to_dggs_array(raster, level, kw...) |> DGGSLayer