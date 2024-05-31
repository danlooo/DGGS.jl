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
    id::Symbol=:layer
)
    level > 0 || error("Level must be positive")

    lon_dim = filter(x -> x isa X || name(x) == lon_name, dims(raster))
    lat_dim = filter(x -> x isa Y || name(x) == lat_name, dims(raster))

    isempty(lon_dim) && error("Longitude dimension not found")
    isempty(lat_dim) && error("Latitude dimension not found")
    lon_dim = lon_dim[1]
    lat_dim = lat_dim[1]

    # data would be upside down if reversed dimensions are used
    issorted(lon_dim) || error("Longitude must be sorted in forward ascending oder")
    issorted(lat_dim) || error("Latitude must be sorted in forward ascending oder")

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

    function raster_to_dggs(xout, xin)
        for (cell_id, cell_indices) in cell_ids_indexlist
            # xout is not a YAXArray anymore
            xout[cell_id.i+1, cell_id.j+1, cell_id.n+1] = agg_func(view(xin, cell_indices))
        end
    end

    cell_array = mapCube(
        raster_to_dggs,
        # mapCube can't find axes of other AbstractDimArrays e.g. Raster
        YAXArray(dims(raster), raster.data),
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

function Makie.plot(a::DGGSArray, args...; type=:globe, kwargs...)
    type == :globe && return plot_globe(a; kwargs...)
    type == :map && return plot_map(a, args...; kwargs...)
    error("Plot type :$type must be one of [:globe, :map]")
end

function plot_globe(a::DGGSArray; resolution::Integer=800)
    longitudes = range(-180, 180; length=resolution * 2)
    latitudes = range(-90, 90; length=resolution)
    geo_array = to_geo_array(a, longitudes, latitudes)

    non_spatial_axes = map(setdiff(DimensionalData.name(geo_array.axes), (:lon, :lat))) do x
        getproperty(geo_array, x)
    end
    min_val = filter_null(minimum)(geo_array)
    max_val = filter_null(maximum)(geo_array)
    min_val == max_val && error("")

    with_theme(theme_black()) do
        fig = Figure()
        ax = fig[1, 1] = LScene(fig[1, 1], show_axis=false)
        side_panel = fig[1, 2] = GridLayout()

        if length(non_spatial_axes) == 0
            # calc texture colors
            texture = map(geo_array |> x -> x[1:length(x.lon), length(x.lat):-1:1]') do val
                try
                    get(ColorSchemes.viridis, val, (min_val, max_val))
                catch
                    RGBA{Float64}(0.15, 0.15, 0.15, 1)
                end
            end
            mesh = Sphere(Point3f(0), 1) |> x -> Tesselation(x, 128) |> uv_mesh
            m_plt = mesh!(
                ax,
                mesh,
                color=texture,
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

            texture = lift(slider_observables...) do slider_values...
                # filter to selected dimensions
                d = Dict()
                for (ax, val) in zip(non_spatial_axes, slider_values)
                    d[name(ax)] = val
                end
                filtered_array = getindex(geo_array; NamedTuple(d)...)

                # calc texture colors
                map(filtered_array |> x -> x[1:length(x.lon), length(x.lat):-1:1]') do val
                    try
                        get(ColorSchemes.viridis, val, (min_val, max_val))
                    catch
                        RGBA{Float64}(0.15, 0.15, 0.15, 1)
                    end
                end
            end
            mesh = Sphere(Point3f(0), 1) |> x -> Tesselation(x, 128) |> uv_mesh
            m_plt = mesh!(
                ax,
                mesh,
                color=texture,
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