function DGGSArray(arr::YAXArray, id::Symbol)
    haskey(arr.properties, "_DGGS") || error("Array is not in DGGS format")

    attrs = arr.properties
    level = attrs["_DGGS"]["level"]
    dggs = DGGSGridSystem(attrs["_DGGS"])

    DGGSArray(arr, attrs, id, level, dggs)
end

function DGGSArray(arr::YAXArray)
    id = get(arr.properties, "name", "layer") |> Symbol
    DGGSArray(arr, id)
end

function DGGSArray(arr::AbstractArray, level::Integer, id=:layer)
    axs = (
        Dim{:q2di_i}(range(0; step=1, length=2^(level - 1))),
        Dim{:q2di_j}(range(0; step=1, length=2^(level - 1))),
        Dim{:q2di_n}(0:11),
    )
    props = Dict{String,Any}(
        "name" => id,
        "_DGGS" => DGGS.Q2DI_DGGS_PROPS
    )
    props["_DGGS"]["level"] = level
    data = YAXArray(axs, arr, props)
    DGGSArray(data)
end

"filter any dimension of a DGGSArray"
Base.getindex(a::DGGSArray, args...; kwargs...) = getindex(a.data, args...; kwargs...) |> DGGSArray

"get a cell of a DGGSArray"
Base.getindex(a::DGGSArray, center::Q2DI; kwargs...) = getindex(a.data, q2di_n=center.n, q2di_i=center.i, q2di_j=center.j; kwargs...)

"get a ring of a DGGArray"
function Base.getindex(a::DGGSArray, center::Q2DI, radius::Integer; kwargs...)
    radius >= 1 || error("radius must not be negative")

    res = a
    if length(kwargs) >= 1
        res = getindex(res; kwargs...)
    end
    res = getindex(res, center, radius, :ring)
    return res
end

"get a disk of a DGGArray"
function Base.getindex(a::DGGSArray, center::Q2DI, range::UnitRange{R}; kwargs...) where {R<:Integer}
    range.start >= 1 || error("Range must start with a positive number")
    range.start == 1 || error("annulus not supported")

    res = a
    if length(kwargs) >= 1
        res = getindex(res; kwargs...)
    end
    res = getindex(res, center, range.stop, :disk)
    return res
end

function Base.getindex(a::DGGSArray, lon::Real, lat::Real, args...; kwargs...)
    center = transform_points([(lon, lat)], a.level)[1]
    res = getindex(a, center, args...; kwargs...)
    return res
end


Base.setindex!(a::DGGSArray, val, i::Q2DI; kwargs...) = Base.setindex!(a.data, val, q2di_n=i.n, q2di_i=i.i, q2di_j=i.j; kwargs...)
Base.setindex!(a::DGGSArray, val, n::Integer, i::Integer, j::Integer; kwargs...) = Base.setindex!(a.data, val, q2di_n=n, q2di_i=i, q2di_j=j; kwargs...)

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

is_spatial(ax) = ax |> name |> String |> x -> startswith(x, "q2di_")

function show_axes(io, axs; hide_if_spatial=true)
    printstyled(io, "Non spatial axes:\n"; color=:white)
    for ax in axs
        hide_if_spatial & is_spatial(ax) && continue

        print(io, "  ")
        printstyled(io, name(ax); color=:red)
        print(io, " $(length(ax)) ")
        print(io, eltype(ax))
        println(io, " points")
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

    println(io, "DGGS:\t\t$(arr.dggs) at level $(arr.level)")
    println(io, "Attributes:\t$(length(arr.attrs))")

    show_axes(io, arr.data.axes)
end

function Base.show(io::IO, arr::DGGSArray)
    "$(arr.id)" |> x -> printstyled(io, x; color=:red)
    print(io, " ")
    if "standard_name" in keys(arr.attrs)
        printstyled(io, arr.attrs["standard_name"]; color=:white)
        print(io, " ")
    end
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

function get_arr_label(a::DGGSArray)
    arr_label = get(a.attrs, "standard_name", "$(a.id)") |> x -> replace(x, "_" => " ")
    if "units" in keys(a.attrs)
        arr_label *= (" (" * a.attrs["units"] * ")")
    end
    return arr_label
end

function to_dggs_array(
    raster::AbstractDimArray,
    level::Integer;
    lon_name::Symbol=:lon,
    lat_name::Symbol=:lat,
    cell_ids::Union{AbstractMatrix,Nothing}=nothing,
    agg_func::Function=filter_null(mean),
    verbose::Bool=true,
    id::Symbol=:layer,
    path::String=tempname()
)
    # Optimized for grids which cell_ids fit in memory

    level > 0 || error("Level must be positive")

    lon_dim = filter(x -> x isa X || name(x) == lon_name, dims(raster))
    lat_dim = filter(x -> x isa Y || name(x) == lat_name, dims(raster))

    isempty(lon_dim) && error("Longitude dimension not found")
    isempty(lat_dim) && error("Latitude dimension not found")
    lon_dim = lon_dim[1]
    lat_dim = lat_dim[1]

    -180 <= minimum(lon_dim) <= maximum(lon_dim) <= 180 || error("$(name(lon_dim)) must be within [-180, 180]")
    -90 <= minimum(lat_dim) <= maximum(lat_dim) <= 90 || error("$(name(lon_dim)) must be within [-90, 90]")

    verbose && @info "Step 1/2: Transform coordinates"
    cell_ids_mat = isnothing(cell_ids) ? transform_points(lon_dim.val, lat_dim.val, level) : cell_ids

    # look up hash map to get list of pixels to be aggregated for each cell
    cell_ids_indexlist = Dict()
    for c in axes(cell_ids_mat, 2)
        for r in axes(cell_ids_mat, 1)
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
        indims=InDims(lon_dim, lat_dim),
        outdims=OutDims(
            Dim{:q2di_i}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_j}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_n}(0:11),
            path=path
        ),
        showprog=true
    ) do xout, xin
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

function to_geo_array(a::DGGSArray, cell_ids::DimArray)
    geo_array = mapCube(
        a.data,
        cell_ids,
        indims=InDims(:q2di_i, :q2di_j, :q2di_n),
        outdims=OutDims(
            Dim{:lon}(dims(cell_ids, :X).val),
            Dim{:lat}(dims(cell_ids, :Y).val)
        )
    ) do xout, xin, cell_ids
        for (i, cell_id) in enumerate(cell_ids)
            xout[i] = xin[cell_id.i+1, cell_id.j+1, cell_id.n+1]
        end
    end

    return YAXArray(geo_array.axes, geo_array.data, a.attrs)
end

function to_geo_array(
    a::DGGSArray,
    longitudes=range(-180, 180; length=800),
    latitudes=range(-90, 90; length=400);
)
    cell_ids = transform_points(longitudes, latitudes, a.level)
    to_geo_array(a, cell_ids)
end

Base.collect(a::DGGSArray) = Base.collect(a.data)

function Makie.plot(a::DGGSArray, args...; type=:globe, kwargs...)
    type == :globe && return plot_globe(a; kwargs...)
    type == :map && return plot_map(a, args...; kwargs...)
    type == :native && return plot_native(a, args...; kwargs...)
    error("Plot type :$type must be one of [:globe, :map, :native]")
end

function plot_globe(a::DGGSArray; resolution::Integer=800)
    # get grid points
    longitudes = range(-180, 180; length=resolution * 2)
    latitudes = range(-90, 90; length=resolution)
    cell_ids = transform_points(longitudes, latitudes, a.level)

    non_spatial_axes = filter(x -> !(name(x) in [:q2di_i, :q2di_j, :q2di_n]), a.data.axes)

    with_theme(theme_black()) do
        fig = Figure()
        ax = fig[1, 1] = LScene(fig[1, 1], show_axis=false)
        side_panel = fig[1, 2] = GridLayout()

        if length(non_spatial_axes) == 0
            geo_array = to_geo_array(a, cell_ids)
            min_val = filter_null(minimum)(geo_array)
            max_val = filter_null(maximum)(geo_array)
            texture = geo_array |>
                      x -> x[1:length(x.lon), length(x.lat):-1:1]' |>
                           collect .|>
                           x -> ismissing(x) ? NaN : x

            msh = Sphere(Point3f(0), 1) |> x -> Tesselation(x, 128) |> uv_mesh
            m_plt = mesh!(
                ax,
                msh,
                color=texture,
                nan_color=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                interpolate=true,
                shading=Makie.NoShading,
            )
            center!(ax.scene)
            cam = Camera3D(
                ax.scene,
                projectiontype=Makie.Perspective,
                cad=true, # prevent dithering
                fixed_axis=true,
                lookat=[0, 0, 0],
                upvector=[0, 0, 1],

                # disable translation
                forward_key=Keyboard.unknown,
                backward_key=Keyboard.unknown,
                up_key=Keyboard.unknown,
                down_key=Keyboard.unknown,
                left_key=Keyboard.unknown,
                right_key=Keyboard.unknown,
                translation_button=Mouse.none,
                pan_left_key=Keyboard.left,
                pan_right_key=Keyboard.right,
                tilt_up_key=Keyboard.up,
                tilt_down_key=Keyboard.down
            )

            north_up_btn = Button(side_panel[2, 1]; label="N↑")
            on(north_up_btn.clicks) do c
                update_cam!(ax.scene, cam.eyeposition[], cam.lookat[], (0.0, 0.0, 1.0))
            end

            # prevent view inside of earth    
            on(cam.eyeposition) do eyeposition
                dist_to_earth_center = norm(cam.eyeposition[])
                min_zoom = 1.1
                if dist_to_earth_center < min_zoom
                    zoom!(ax.scene, min_zoom)
                end
            end

            # north_up_btn.labelcolor = :red
            north_up_btn.buttoncolor = :black
            north_up_btn.buttoncolor_hover = :black
            north_up_btn.labelcolor_hover = :white

            cb = Colorbar(side_panel[1, 1]; limits=(min_val, max_val), label=get_arr_label(a))
            fig
        else
            sliders = map(non_spatial_axes) do ax
                (
                    label=ax |> name |> String,
                    range=1:length(ax),
                    format=x -> "$(ax[x])",
                    color_inactive=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                    color_active=:white,
                    color_active_dimmed=:white
                )
            end
            slider_grid = SliderGrid(fig[2, 1], sliders...)
            slider_observables = [s.value for s in slider_grid.sliders]

            geo_array = lift(slider_observables...) do slider_values...
                # filter to selected dimensions
                d = Dict()
                for (ax, val) in zip(non_spatial_axes, slider_values)
                    d[name(ax)] = val
                end

                getindex(a; NamedTuple(d)...) |> x -> to_geo_array(x, cell_ids)
            end

            min_val = @lift filter_null(minimum)($geo_array)
            max_val = @lift filter_null(maximum)($geo_array)

            texture = @lift $geo_array |>
                            x -> x[1:length(x.lon), length(x.lat):-1:1]' |>
                                 collect .|>
                                 x -> ismissing(x) ? NaN : x

            msh = Sphere(Point3f(0), 1) |> x -> Tesselation(x, 128) |> uv_mesh
            m_plt = mesh!(
                ax,
                msh,
                color=texture,
                nan_color=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                interpolate=true,
                shading=Makie.NoShading,
            )
            center!(ax.scene)
            cam = Camera3D(
                ax.scene,
                projectiontype=Makie.Perspective,
                cad=true, # prevent dithering
                fixed_axis=true,
                lookat=[0, 0, 0],
                upvector=[0, 0, 1],

                # disable translation
                forward_key=Keyboard.unknown,
                backward_key=Keyboard.unknown,
                up_key=Keyboard.unknown,
                down_key=Keyboard.unknown,
                left_key=Keyboard.unknown,
                right_key=Keyboard.unknown,
                translation_button=Mouse.none,
                pan_left_key=Keyboard.left,
                pan_right_key=Keyboard.right,
                tilt_up_key=Keyboard.up,
                tilt_down_key=Keyboard.down
            )

            north_up_btn = Button(side_panel[2, 1]; label="N↑")
            on(north_up_btn.clicks) do c
                update_cam!(ax.scene, cam.eyeposition[], cam.lookat[], (0.0, 0.0, 1.0))
            end

            # prevent view inside of earth    
            on(cam.eyeposition) do eyeposition
                dist_to_earth_center = norm(cam.eyeposition[])
                min_zoom = 1.1
                if dist_to_earth_center < min_zoom
                    zoom!(ax.scene, min_zoom)
                end
            end

            # north_up_btn.labelcolor = :red
            north_up_btn.buttoncolor = :black
            north_up_btn.buttoncolor_hover = :black
            north_up_btn.labelcolor_hover = :white

            cb = Colorbar(side_panel[1, 1]; limits=@lift(($min_val, $max_val)), label=get_arr_label(a))
            fig
        end
    end
end

function plot_map(
    a::DGGSArray;
    longitudes=range(-180, 180; length=800),
    latitudes=range(-90, 90; length=400)
)
    #TODO: convert to geo_array lazyly only after slider values were selected
    geo_array = to_geo_array(a, longitudes, latitudes)
    non_spatial_axes = map(setdiff(DimensionalData.name(geo_array.axes), (:lon, :lat))) do x
        getproperty(geo_array, x)
    end

    with_theme(theme_black()) do
        fig = Figure()
        heatmap_ax = (
            backgroundcolor=RGBA{Float64}(0.15, 0.15, 0.15, 1),
            aspect=1
        )

        if length(non_spatial_axes) == 0
            h, heatmap_plt = heatmap(fig[1, 1], geo_array, axis=heatmap_ax)
            cb = Colorbar(fig[1, 2], heatmap_plt; label=get_arr_label(a))
            fig
        else
            sliders = map(non_spatial_axes) do ax
                (
                    label=ax |> name |> String,
                    range=1:length(ax),
                    format=x -> "$(ax[x])",
                    color_inactive=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                    color_active=:white,
                    color_active_dimmed=:white
                )
            end
            slider_grid = SliderGrid(fig[2, 1], sliders...)
            slider_observables = [s.value for s in slider_grid.sliders]

            texture = lift(slider_observables...) do slider_values...
                d = Dict()
                for (ax, val) in zip(non_spatial_axes, slider_values)
                    d[name(ax)] = val
                end
                getindex(geo_array; NamedTuple(d)...)
            end

            h, heatmap_plt = heatmap(fig[1, 1], texture, axis=heatmap_ax)
            cb = Colorbar(fig[1, 2], heatmap_plt; label=get_arr_label(a))
            fig
        end
    end
end

"""
Plot a DGGSArray nativeley on a icosahedron
"""
function plot_native(a::DGGSArray)
    non_spatial_axes = filter(x -> !(name(x) in [:q2di_i, :q2di_j, :q2di_n]), a.data.axes)

    with_theme(theme_black()) do
        fig = Figure()
        ax = fig[1, 1] = LScene(fig[1, 1], show_axis=false)
        side_panel = fig[1, 2] = GridLayout()

        if length(non_spatial_axes) == 0
            native_array = a.data[q2di_n=2:11] |> collect # ignore 2 vertices at quad 1 and 12
            min_val = filter_null(minimum)(native_array)
            max_val = filter_null(maximum)(native_array)
            texture = native_array |>
                      x -> reshape(x, size(x)[1], 10 * size(x)[1]) .|>
                           x -> ismissing(x) ? NaN : x
            msh = load(artifact"isea-obj" * "/isea.obj")
            m_plt = mesh!(
                ax,
                msh,
                color=texture,
                nan_color=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                interpolate=true,
                shading=Makie.NoShading,
            )
            center!(ax.scene)
            cam = Camera3D(
                ax.scene,
                projectiontype=Makie.Perspective,
                cad=true, # prevent dithering
                fixed_axis=true,
                lookat=[0, 0, 0],
                upvector=[0, 0, 1],

                # disable translation
                forward_key=Keyboard.unknown,
                backward_key=Keyboard.unknown,
                up_key=Keyboard.unknown,
                down_key=Keyboard.unknown,
                left_key=Keyboard.unknown,
                right_key=Keyboard.unknown,
                translation_button=Mouse.none,
                pan_left_key=Keyboard.left,
                pan_right_key=Keyboard.right,
                tilt_up_key=Keyboard.up,
                tilt_down_key=Keyboard.down
            )

            north_up_btn = Button(side_panel[2, 1]; label="N↑")
            on(north_up_btn.clicks) do c
                update_cam!(ax.scene, cam.eyeposition[], cam.lookat[], (0.0, 0.0, 1.0))
            end

            # prevent view inside of earth    
            on(cam.eyeposition) do eyeposition
                dist_to_earth_center = norm(cam.eyeposition[])
                min_zoom = 1.1
                if dist_to_earth_center < min_zoom
                    zoom!(ax.scene, min_zoom)
                end
            end

            # north_up_btn.labelcolor = :red
            north_up_btn.buttoncolor = :black
            north_up_btn.buttoncolor_hover = :black
            north_up_btn.labelcolor_hover = :white

            cb = Colorbar(side_panel[1, 1]; limits=(min_val, max_val), label=get_arr_label(a))
            fig
        else
            sliders = map(non_spatial_axes) do ax
                (
                    label=ax |> name |> String,
                    range=1:length(ax),
                    format=x -> "$(ax[x])",
                    color_inactive=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                    color_active=:white,
                    color_active_dimmed=:white
                )
            end
            slider_grid = SliderGrid(fig[2, 1], sliders...)
            slider_observables = [s.value for s in slider_grid.sliders]

            native_array = lift(slider_observables...) do slider_values...
                # filter to selected dimensions
                d = Dict()
                for (ax, val) in zip(non_spatial_axes, slider_values)
                    d[name(ax)] = val
                end

                getindex(a; NamedTuple(d)...) |> x -> x.data[q2di_n=2:11] |> collect
            end

            min_val = @lift filter_null(minimum)($native_array)
            max_val = @lift filter_null(maximum)($native_array)

            texture = @lift $native_array |>
                            x -> reshape(x, size(x)[1], 10 * size(x)[1]) .|>
                                 x -> ismissing(x) ? NaN : x

            msh = load(artifact"isea-obj" * "/isea.obj")
            m_plt = mesh!(
                ax,
                msh,
                color=texture,
                nan_color=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                interpolate=true,
                shading=Makie.NoShading,
            )
            center!(ax.scene)
            cam = Camera3D(
                ax.scene,
                projectiontype=Makie.Perspective,
                cad=true, # prevent dithering
                fixed_axis=true,
                lookat=[0, 0, 0],
                upvector=[0, 0, 1],

                # disable translation
                forward_key=Keyboard.unknown,
                backward_key=Keyboard.unknown,
                up_key=Keyboard.unknown,
                down_key=Keyboard.unknown,
                left_key=Keyboard.unknown,
                right_key=Keyboard.unknown,
                translation_button=Mouse.none,
                pan_left_key=Keyboard.left,
                pan_right_key=Keyboard.right,
                tilt_up_key=Keyboard.up,
                tilt_down_key=Keyboard.down
            )

            north_up_btn = Button(side_panel[2, 1]; label="N↑")
            on(north_up_btn.clicks) do c
                update_cam!(ax.scene, cam.eyeposition[], cam.lookat[], (0.0, 0.0, 1.0))
            end

            # prevent view inside of earth    
            on(cam.eyeposition) do eyeposition
                dist_to_earth_center = norm(cam.eyeposition[])
                min_zoom = 1.1
                if dist_to_earth_center < min_zoom
                    zoom!(ax.scene, min_zoom)
                end
            end

            # north_up_btn.labelcolor = :red
            north_up_btn.buttoncolor = :black
            north_up_btn.buttoncolor_hover = :black
            north_up_btn.labelcolor_hover = :white

            cb = Colorbar(side_panel[1, 1]; limits=@lift(($min_val, $max_val)), label=get_arr_label(a))
            fig
        end
    end

end

#
# Arithmetics
#

import Base: broadcasted, +, -, *, /, \

function Base.broadcasted(f::Function, a::DGGSArray, number::Real)
    data = f.(a.data, number)

    # meta data may be invalidated after transformation
    properties = Dict("_DGGS" => deepcopy(data.properties["_DGGS"]))
    array = YAXArray(data.axes, data.data, properties)

    res = DGGSArray(array)
    return res
end

for f in (:/, :\, :*)
    if f !== :/
        @eval ($f)(number::Number, a::DGGSArray) = Base.broadcasted($f, a, number)
    end
    if f !== :\
        @eval ($f)(a::DGGSArray, number::Number) = Base.broadcasted($f, a, number)
    end
end
