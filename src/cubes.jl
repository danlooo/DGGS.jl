function Base.show(io::IO, ::MIME"text/plain", cube::DGGSArray)
    println(io, "DGGSArray at level $(cube.level)")
    Base.show(io, "text/plain", cube.data.axes)
end

function Base.getindex(cell_cube::DGGSArray, i::Q2DI)
    cell_cube.data[q2di_n=At(i.n), q2di_i=At(i.i), q2di_j=At(i.j)]
end

function Base.getindex(cell_cube::DGGSArray, lon::Real, lat::Real)
    cell_id = _transform_points(lon, lat, cell_cube.level)[1, 1]
    cell_cube.data[q2di_n=At(cell_id.n), q2di_i=At(cell_id.i), q2di_j=At(cell_id.j)]
end

"Apply function f after filtering of missing and NAN values"
function filter_null(f)
    x -> x |> filter(!ismissing) |> filter(!isnan) |> f
end

function map_geo_to_cell_cube(xout, xin, cell_ids_indexlist, agg_func)
    for (cell_id, cell_indices) in cell_ids_indexlist
        # xout is not a YAXArray anymore
        xout[cell_id.i+1, cell_id.j+1, cell_id.n+1] = agg_func(view(xin, cell_indices))
    end
end


"maximial i or j value in Q2DI index given a level"
max_ij(level) = level <= 3 ? level - 1 : 2^(level - 2)

function to_dggs_array(raster::AbstractDimArray, level::Integer; agg_func::Function=filter_null(mean), cell_ids::Union{AbstractMatrix,Nothing}=nothing, lon_name::Symbol=:lon, lat_name::Symbol=:lat)
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

    @info "Step 1/2: Transform coordinates"
    cell_ids_mat = isnothing(cell_ids) ? transform_points(lon_dim.val, lat_dim.val, level) : cell_ids

    # TODO: what if e.g. lon goes from 0 to 365?
    # TODO: Reverse inverted geo axes

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

    @info "Step 2/2: Re-grid the data"
    cell_cube = mapCube(
        map_geo_to_cell_cube,
        # mapCube can't find axes of other AbstractDimArrays e.g. Raster
        YAXArray(dims(raster), raster.data),
        cell_ids_indexlist,
        agg_func,
        indims=InDims(:X, :Y),
        outdims=OutDims(
            Dim{:q2di_i}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_j}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_n}(0:11)
        ),
        showprog=true
    )
    cell_cube = YAXArray(cell_cube.axes, cell_cube.data, metadata(raster))
    DGGSArray(cell_cube, level)
end

function to_dggs_array(raster::AbstractMatrix, lon_range::AbstractVector, lat_range::AbstractVector, level::Integer; kwargs...)
    raster = DimArray(raster, (X(lon_range), Y(lat_range)))
    return to_dggs_array(raster, level; kwargs...)
end

Base.getindex(cell_cube::DGGSArray; i...) = Base.getindex(cell_cube.data; i...) |> x -> DGGSArray(x, cell_cube.level)

function map_cell_to_geo_cube(xout, xin, cell_ids_mat)
    for (i, cell_id) in enumerate(cell_ids_mat)
        xout[i] = xin[cell_id.i+1, cell_id.j+1, cell_id.n+1]
    end
end

function to_geo_cube(cell_cube::DGGSArray; longitudes=-180:180, latitudes=-90:90)
    cell_ids_mat = transform_points(longitudes, latitudes, cell_cube.level).data

    geo_array = mapCube(
        map_cell_to_geo_cube,
        cell_cube.data,
        cell_ids_mat,
        indims=InDims(:q2di_i, :q2di_j, :q2di_n),
        outdims=OutDims(
            Dim{:lon}(longitudes),
            Dim{:lat}(latitudes)
        )
    )
    return geo_array
end


function to_geo_cube(cell_cube::DGGSArray, cell_ids::DimArray{Q2DI{T},2}) where {T<:Integer}
    geo_array = mapCube(
        map_cell_to_geo_cube,
        cell_cube.data,
        cell_ids.data,
        indims=InDims(:q2di_i, :q2di_j, :q2di_n),
        outdims=OutDims(dims(cell_ids)...)
    )
    return geo_array
end


function color_value(value, color_scale::ColorScale; null_color=RGBA{Float64}(0.15, 0.15, 0.15, 1))
    ismissing(value) && return null_color
    isnan(value) && return null_color
    return color_scale.schema[value] |> RGBA
end

function get_non_spatial_cube_axes(cell_cube)
    non_spatial_cube_axes = []
    for (i, ax) in enumerate(cell_cube.data.axes)
        name(ax) in [:q2di_i, :q2di_j, :q2di_n] && continue

        entry = Dict(
            :slider => (
                label=ax |> name |> String,
                range=1:length(ax),
                format=x -> "$(ax[x])",
                color_inactive=RGBA{Float64}(0.15, 0.15, 0.15, 1),
                color_active=:white,
                color_active_dimmed=:white
            ),
            :dim => ax
        )
        push!(non_spatial_cube_axes, entry)
    end
    non_spatial_cube_axes
end

function plot_geo(cell_cube::DGGSArray; resolution::Real=800)
    longitudes = range(-180, 180, length=resolution)
    latitudes = range(-90, 90, length=resolution)
    cell_ids = transform_points(longitudes, latitudes, cell_cube.level)
    plot_geo(cell_cube, cell_ids)
end

function plot_geo(cell_cube::DGGSArray, cell_ids::DimArray{Q2DI{T},2}) where {T<:Integer}
    longitudes = dims(cell_ids, :X)
    latitudes = dims(cell_ids, :Y)

    with_theme(theme_black()) do
        fig = Figure()
        ax = fig[1, 1] = LScene(fig[1, 1], show_axis=false)
        side_panel = fig[1, 2] = GridLayout()
        cb = Colorbar(side_panel[1, 1])

        non_spatial_cube_axes = get_non_spatial_cube_axes(cell_cube)

        if length(non_spatial_cube_axes) > 0
            slider_grid = SliderGrid(fig[2, 1], [x[:slider] for x in non_spatial_cube_axes]...)
            slider_observables = [s.value for s in slider_grid.sliders]

            texture = lift(slider_observables...) do slider_values...
                d = Dict()
                for (ax, val) in zip(non_spatial_cube_axes, slider_values)
                    d[name(ax[:dim])] = val
                end
                d = NamedTuple(d)
                filtered_cell_cube = getindex(cell_cube; d...)

                geo_cube = to_geo_cube(filtered_cell_cube, cell_ids)
                color_scale = ColorScale(ColorSchemes.viridis, filter_null(minimum)(geo_cube) |> floor |> Int, filter_null(maximum)(geo_cube) |> ceil |> Int)
                cb.limits[] = (color_scale.min_value, color_scale.max_value)
                texture = map(x -> color_value(x, color_scale), geo_cube[1:length(longitudes), length(latitudes):-1:1]')
            end
            texture
        else
            geo_cube = to_geo_cube(cell_cube, cell_ids)
            color_scale = ColorScale(ColorSchemes.viridis, filter_null(minimum)(geo_cube.data), filter_null(maximum)(geo_cube.data))
            texture = map(x -> color_value(x, color_scale), geo_cube.data[1:length(longitudes), length(latitudes):-1:1]')
            texture
        end

        mesh = Sphere(Point3f(0), 1) |> x -> Tesselation(x, 128) |> uv_mesh
        mesh!(
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
        fig
    end
end

function plot_geo(cell_cube::DGGSArray, bbox::HyperRectangle{2,Float32}; resolution::Int64=800)
    min_x, min_y = bbox.origin
    max_x = min_x + bbox.widths[1]
    max_y = min_y + bbox.widths[2]
    longitudes = range(min_x, max_x; length=resolution)
    latitudes = range(min_y, max_y; length=resolution)
    cell_ids = transform_points(longitudes, latitudes, cell_cube.level)
    geo_cube = to_geo_cube(cell_cube, cell_ids)


    with_theme(theme_black()) do
        fig = Figure()
        non_spatial_cube_axes = get_non_spatial_cube_axes(cell_cube)

        if length(non_spatial_cube_axes) > 0
            slider_grid = SliderGrid(fig[2, 1], [x[:slider] for x in non_spatial_cube_axes]...)
            slider_observables = [s.value for s in slider_grid.sliders]

            geo_cube = lift(slider_observables...) do slider_values...
                d = Dict()
                for (ax, val) in zip(non_spatial_cube_axes, slider_values)
                    d[name(ax[:dim])] = val
                end
                filtered_cell_cube = getindex(cell_cube; NamedTuple(d)...)
                geo_cube = to_geo_cube(filtered_cell_cube, cell_ids)
                geo_cube
            end
            geo_cube
        else
            geo_cube = Observable(to_geo_cube(cell_cube, cell_ids))
        end

        min_value = filter_null(minimum)(geo_cube[].data) |> floor |> Int
        max_value = filter_null(maximum)(geo_cube[].data) |> ceil |> Int
        color_scale = ColorScale(ColorSchemes.viridis, min_value, max_value)
        cb = Colorbar(fig[1, 2])
        cb.limits[] = (color_scale.min_value, color_scale.max_value)

        texture = @lift $geo_cube |> DimArray
        heatmap(fig[1, 1], texture, axis=(backgroundcolor=RGBA{Float64}(0.15, 0.15, 0.15, 1), aspect=1))
        fig
    end
end

function plot_native(cell_cube::DGGSArray)
    cell_cube = cell_cube[q2di_n=2:11] # ignore 2 vertices at quad 1 and 12

    with_theme(theme_black()) do
        fig = Figure()
        ax = fig[1, 1] = LScene(fig[1, 1], show_axis=false)
        side_panel = fig[1, 2] = GridLayout()
        cb = Colorbar(side_panel[1, 1])
        non_spatial_cube_axes = get_non_spatial_cube_axes(cell_cube)

        if length(non_spatial_cube_axes) > 0
            slider_grid = SliderGrid(fig[2, 1], [x[:slider] for x in non_spatial_cube_axes]...)
            slider_observables = [s.value for s in slider_grid.sliders]

            texture = lift(slider_observables...) do slider_values...
                d = Dict()
                for (ax, val) in zip(non_spatial_cube_axes, slider_values)
                    d[name(ax[:dim])] = val
                end
                d = NamedTuple(d)
                filtered_cell_cube = getindex(cell_cube; d...)

                color_scale = ColorScale(ColorSchemes.viridis, filter_null(minimum)(filtered_cell_cube.data) |> floor |> Int, filter_null(maximum)(filtered_cell_cube.data) |> ceil |> Int)
                cb.limits[] = (color_scale.min_value, color_scale.max_value)

                texture = filtered_cell_cube.data
                texture = Array(texture)
                texture = reshape(texture, size(texture)[1], 10 * size(texture)[1])
                texture = map(x -> color_value(x, color_scale), texture)
                texture
            end
        else
            color_scale = ColorScale(ColorSchemes.viridis, filter_null(minimum)(cell_cube.data) |> floor |> Int, filter_null(maximum)(cell_cube.data) |> ceil |> Int)
            texture = Array(cell_cube.data)
            texture = reshape(texture, size(texture)[1], 10 * size(texture)[1])
            texture = map(x -> color_value(x, color_scale), texture)
            texture
        end

        cam = Camera3D(
            ax.scene,
            show_axis=false,
            projectiontype=Makie.Perspective,
            clipping_mode=:static,
            cad=true, # prevent dithering
            lookat=[0, 0, 0],
            upvector=[0, 0, 1],
            fixed_axis=true,
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
        msh = load(artifact"isea-obj" * "/isea.obj")
        mesh!(ax, msh, color=texture, shading=NoShading)
        fig
    end
end

function Makie.plot(cell_cube::DGGSArray, args...; type=:geo, kwargs...)
    if type == :geo
        plot_geo(cell_cube, args...; kwargs...)
    elseif type == :native
        plot_native(cell_cube, args...; kwargs...)
    end
end