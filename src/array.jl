function DGGSArray(arr::YAXArray, id=:layer)
    haskey(arr.properties, "_DGGS") || error("Array is not in DGGS format")

    attrs = arr.properties
    level = attrs["_DGGS"]["level"]
    dggs = DGGSGridSystem(attrs["_DGGS"])

    DGGSArray(arr, attrs, id, level, dggs)
end

function show_nonspatial_axes(io::IO, arr::DGGSArray)
    non_spatial_axes = filter(x -> !startswith(String(x), "q2di"), DimensionalData.name(arr.data.axes))
    if length(non_spatial_axes) == 1
        printstyled(io, "(:"; color=:white)
        printstyled(io, non_spatial_axes[1]; color=:white)
        printstyled(io, ") "; color=:white)
    elseif length(non_spatial_axes) > 1
        printstyled(io, non_spatial_axes; color=:white)
        printstyled(io, " "; color=:white)
    end
end

function Base.show(io::IO, ::MIME"text/plain", arr::DGGSArray)
    printstyled(io, typeof(arr); color=:white)
    println(io, "")
    println(io, "Name:\t\t$(arr.id)")

    if "title" in keys(arr.attrs)
        println(io, "Title:\t\t$(arr.attrs["title"])")
    end

    if "standard_name" in keys(arr.attrs)
        println(io, "Standard name:\t$(arr.attrs["standard_name"])")
    end
    if "units" in keys(arr.attrs)
        println(io, "Units:\t\t$(arr.attrs["units"])")
    end

    println(io, "DGGS:\t\t$(arr.dggs)")
    println(io, "Level:\t\t$(arr.level)")
    println(io, "Attributes: $(length(arr.attrs))")

    println(io, "Non spatial axes:")
    for ax in arr.data.axes
        ax_name = DimensionalData.name(ax)
        startswith(String(ax_name), "q2di") && continue

        print(io, "  ")
        printstyled(io, ax_name; color=:red)
        print(io, " ")
        print(io, eltype(ax))
        println(io, "")
    end
end

function Base.show(io::IO, arr::DGGSArray)
    "$(arr.id)" |> x -> printstyled(io, x; color=:red)
    print(io, " ")
    get(arr.attrs, "standard_name", "") |> x -> printstyled(io, x; color=:white)
    print(io, " ")
    show_nonspatial_axes(io, arr)
    get(arr.attrs, "units", "") |> x -> printstyled(io, x; color=:blue)
    print(io, " ")
    eltype(arr.data) |> x -> print(io, x)
    print(io, " ")
    "cell_methods" in keys(arr.attrs) && print(io, "aggregated")
end

"Apply function f after filtering of missing and NAN values"
function filter_null(f)
    x -> x |> filter(!ismissing) |> filter(!isnan) |> f
end

function to_dggs_array(
    raster::AbstractDimArray,
    level::Integer;
    lon_name::Symbol=:lon,
    lat_name::Symbol=:lat,
    cell_ids::Union{AbstractMatrix,Nothing}=nothing,
    agg_func::Function=filter_null(mean),
    verbose::Bool=true,
    id::Symbol=:layer
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
    ) do xout, xin, cell_ids_indexlist, agg_func
        for (cell_id, cell_indices) in cell_ids_indexlist
            # xout is not a YAXArray anymore
            xout[cell_id.i+1, cell_id.j+1, cell_id.n+1] = agg_func(view(xin, cell_indices))
        end
    end

    props = raster |> metadata |> typeof == Dict{String,Any} ? deepcopy(metadata(raster)) : Dict{String,Any}()
    props["_DGGS"] = deepcopy(Q2DI_DGGS_PROPS)
    props["_DGGS"]["level"] = level
    cell_array = YAXArray(cell_array.axes, cell_array.data, props)
    DGGSArray(cell_array, id)
end

function to_geo_array(
    a::DGGSArray,
    longitudes=range(-180, 180; length=800),
    latitudes=range(-90, 90; length=400);
    cell_ids=nothing
)
    if isnothing(cell_ids)
        cell_ids = transform_points(longitudes, latitudes, a.level)
    end

    geo_array = mapCube(
        a.data,
        cell_ids,
        indims=InDims(:q2di_i, :q2di_j, :q2di_n),
        outdims=OutDims(
            Dim{:lon}(longitudes),
            Dim{:lat}(latitudes)
        )
    ) do xout, xin, cell_ids
        for (i, cell_id) in enumerate(cell_ids)
            xout[i] = xin[cell_id.i+1, cell_id.j+1, cell_id.n+1]
        end
    end

    return YAXArray(geo_array.axes, geo_array.data, a.attrs)
end

function Makie.plot(
    a::DGGSArray;
    longitudes=range(-180, 180; length=800),
    latitudes=range(-90, 90; length=400)
)
    geo_array = to_geo_array(a, longitudes, latitudes)
    non_spatial_axes = setdiff(DimensionalData.name(geo_array.axes), (:lon, :lat))

    fig = Figure()
    if length(non_spatial_axes) == 0
        with_theme(theme_black()) do
            heatmap(fig[1, 1], geo_array, axis=(backgroundcolor=RGBA{Float64}(0.15, 0.15, 0.15, 1), aspect=1))
            fig
        end
    else
        sliders = [(
            label=String(key),
            range=1:length(getproperty(geo_array, key))
        ) for key in non_spatial_axes]
        slider_grid = SliderGrid(fig[2, 1], sliders...)

        with_theme(theme_black()) do
            fig = Figure()
            heatmap(fig[1, 1], geo_array, axis=(backgroundcolor=RGBA{Float64}(0.15, 0.15, 0.15, 1), aspect=1))
            fig
        end
    end
end