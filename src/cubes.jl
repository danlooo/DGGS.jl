function Base.show(io::IO, ::MIME"text/plain", cube::CellCube)
    println(io, "DGGS CellCube at level $(cube.level)")
    Base.show(io, "text/plain", cube.data.axes)
end

function CellCube(path::String, lon_dim, lat_dim, level)
    geo_cube = GeoCube(path::String, lon_dim, lat_dim)
    CellCube(geo_cube, level)
end

function GeoCube(path::String, lon_dim, lat_dim)
    array = Cube(path)
    array = renameaxis!(array, lon_dim => :lon)
    array = renameaxis!(array, lat_dim => :lat)

    -180 <= minimum(array.lon) < maximum(array.lon) <= 180 || error("Longitudes must be within [-180, 180]")
    -90 <= minimum(array.lat) < maximum(array.lat) <= 90 || error("Longitudes must be within [-180, 180]")

    GeoCube(array)
end

function GeoCube(data::AbstractMatrix{<:Number}, lon_range::AbstractRange{<:Real}, lat_range::AbstractRange{<:Real})
    axlist = (
        Dim{:lat}(lat_range),
        Dim{:lon}(lon_range)
    )
    geo_array = YAXArray(axlist, data)
    geo_cube = GeoCube(geo_array)
    return geo_cube
end

function Base.show(io::IO, ::MIME"text/plain", cube::GeoCube)
    println(io, "DGGS GeoCube")
    Base.show(io, "text/plain", cube.data.axes)
end

function Base.getindex(cell_cube::CellCube, i::Q2DI)
    cell_cube.data[q2di_n=At(i.n), q2di_i=At(i.i), q2di_j=At(i.j)]
end

function Base.getindex(cell_cube::CellCube, lon::Real, lat::Real)
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

function CellCube(path::String, level; kwargs...)
    geo_cube = GeoCube(path)
    cell_cube = CellCube(geo_cube, level; kwargs...)
    return cell_cube
end

"maximial i or j value in Q2DI index given a level"
max_ij(level) = level <= 3 ? level - 1 : 2^(level - 2)

function CellCube(geo_cube::GeoCube, level, agg_func=filter_null(mean))
    @info "Step 1/2: Transform coordinates"
    cell_ids_mat = transform_points(geo_cube.data.lon, geo_cube.data.lat, level)

    cell_ids_indexlist = Dict()
    for c in 1:size(cell_ids_mat, 2)
        for r in 1:size(cell_ids_mat, 1)
            cell_id = cell_ids_mat[r, c]
            if cell_id in keys(cell_ids_indexlist)
                push!(cell_ids_indexlist[cell_id], CartesianIndex(r, c))
            else
                cell_ids_indexlist[cell_id] = [CartesianIndex(r, c)]
            end
        end
    end

    @info "Step 2/2: Re-grid the data"
    cell_cube = mapCube(
        map_geo_to_cell_cube,
        geo_cube.data,
        cell_ids_indexlist,
        agg_func,
        indims=InDims(:lon, :lat),
        outdims=OutDims(
            Dim{:q2di_i}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_j}(range(0; step=1, length=2^(level - 1))),
            Dim{:q2di_n}(0:11)
        ),
        showprog=true
    )
    return CellCube(cell_cube, level)
end

Base.getindex(cell_cube::CellCube; i...) = Base.getindex(cell_cube.data; i...) |> x -> CellCube(x, cell_cube.level)

function map_cell_to_geo_cube(xout, xin, cell_ids_mat)
    for (i, cell_id) in enumerate(cell_ids_mat)
        xout[i] = xin[cell_id.i+1, cell_id.j+1, cell_id.n+1]
    end
end

function GeoCube(cell_cube::CellCube; longitudes=-180:180, latitudes=-90:90, cell_ids_mat=nothing)
    # transforming points is the slowest step
    # re-load from cache or re-use for all time points
    if isnothing(cell_ids_mat)
        cell_ids_mat = transform_points(longitudes, latitudes, cell_cube.level)
    end

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
    return GeoCube(geo_array)
end

function color_value(value, color_scale::ColorScale; null_color=RGBA{Float64}(0.15, 0.15, 0.15, 1))
    ismissing(value) && return null_color
    isnan(value) && return null_color
    return color_scale.schema[value] |> RGBA
end

function Makie.plot(cell_cube::CellCube; resolution::Int64=800)
    # texture for plot in equirectangular geographic lat/lon projection
    longitudes = range(-180, 180, length=resolution * 2)
    latitudes = range(-90, 90, length=resolution)
    cell_ids_mat = transform_points(longitudes, latitudes, cell_cube.level)

    with_theme(theme_black()) do
        fig = Figure()
        ax = fig[1, 1] = LScene(fig[1, 1], show_axis=false)
        side_panel = fig[1, 2] = GridLayout()
        cb = Colorbar(side_panel[1, 1])

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

                geo_cube = GeoCube(filtered_cell_cube; longitudes, latitudes, cell_ids_mat)
                color_scale = ColorScale(ColorSchemes.viridis, filter_null(minimum)(geo_cube.data.data) |> floor |> Int, filter_null(maximum)(geo_cube.data.data) |> ceil |> Int)
                cb.limits[] = (color_scale.min_value, color_scale.max_value)
                texture = map(x -> color_value(x, color_scale), geo_cube.data.data[1:length(longitudes), length(latitudes):-1:1]')
            end
            texture
        else
            geo_cube = GeoCube(cell_cube; longitudes, latitudes, cell_ids_mat)
            color_scale = ColorScale(ColorSchemes.viridis, filter_null(minimum)(geo_cube.data.data), filter_null(maximum)(geo_cube.data.data))
            texture = map(x -> color_value(x, color_scale), geo_cube.data.data[1:length(longitudes), length(latitudes):-1:1]')
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

function Makie.plot(cell_cube::CellCube, bbox::HyperRectangle{2,Float32}; resolution::Int64=800)
    min_x, min_y = bbox.origin
    max_x = min_x + bbox.widths[1]
    max_y = min_y + bbox.widths[2]
    longitudes = range(min_x, max_x; length=resolution)
    latitudes = range(min_y, max_y; length=resolution)
    cell_ids_mat = transform_points(longitudes, latitudes, cell_cube.level)
    geo_cube = GeoCube(cell_cube; longitudes, latitudes, cell_ids_mat)


    with_theme(theme_black()) do
        fig = Figure()
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

        if length(non_spatial_cube_axes) > 0
            slider_grid = SliderGrid(fig[2, 1], [x[:slider] for x in non_spatial_cube_axes]...)
            slider_observables = [s.value for s in slider_grid.sliders]

            geo_cube = lift(slider_observables...) do slider_values...
                d = Dict()
                for (ax, val) in zip(non_spatial_cube_axes, slider_values)
                    d[name(ax[:dim])] = val
                end
                filtered_cell_cube = getindex(cell_cube; NamedTuple(d)...)
                geo_cube = GeoCube(filtered_cell_cube; longitudes, latitudes, cell_ids_mat)
                geo_cube
            end
            geo_cube
        else
            geo_cube = Observable(GeoCube(cell_cube; longitudes, latitudes, cell_ids_mat))
        end

        min_value = filter_null(minimum)(geo_cube[].data.data) |> floor |> Int
        max_value = filter_null(maximum)(geo_cube[].data.data) |> ceil |> Int
        color_scale = ColorScale(ColorSchemes.viridis, min_value, max_value)
        cb = Colorbar(fig[1, 2])
        cb.limits[] = (color_scale.min_value, color_scale.max_value)

        texture = @lift $geo_cube.data |> DimArray
        heatmap(fig[1, 1], texture, axis=(backgroundcolor=RGBA{Float64}(0.15, 0.15, 0.15, 1), aspect=1))
        fig
    end
end