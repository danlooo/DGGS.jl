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

"""
Filter CellCube using a query string

Name and value of the dimension to be filtered are seperated by `=`.
Multiple filters may be provided. They must be separated by `,`.
All elements of the resulting CellCube passed all given filters.

Examples: `all`, `time=2023-11-04T00:00:00,band=4`
"""
function query(cell_cube::CellCube, query_str::String="all")
    query_str in ["all", ""] && return cell_cube

    query_d = Dict()
    # Try parsing these types with descending priority
    prefered_types = [DateTime, Int64, Float64, Bool, String]

    for dim_query in split(query_str, ",")
        k, v = split(dim_query, "=")

        for t in prefered_types
            v = tryparse(t, v)
            if isnothing(v)
                continue
            end
            query_d[Symbol(k)] = At(v)
            break
        end
    end

    data = Base.getindex(cell_cube.data; NamedTuple(query_d)...)
    return CellCube(data, cell_cube.level)
end

function Base.mapslices(agg_func, cell_cube::CellCube; dims, kw...)
    data = mapslices(agg_func, cell_cube.data; dims, kw...)
    return CellCube(data, cell_cube.level)
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

function plot(cell_cube::CellCube; resolution::Int64=800)
    cell_cube.data.axes .|> name == (:q2di_i, :q2di_j, :q2di_n) || error("cell_cube has Too many dimensions for plotting. Please consider function query to filter.")

    longitudes = range(-180, 180, length=resolution * 2)
    latitudes = range(-90, 90, length=resolution)

    geo_cube = GeoCube(cell_cube; longitudes, latitudes)
    color_scale = ColorScale(ColorSchemes.viridis, filter_null(minimum)(geo_cube.data.data), filter_null(maximum)(geo_cube.data.data))
    texture = map(x -> color_value(x, color_scale), geo_cube.data.data[1:length(longitudes), length(latitudes):-1:1]')
    mesh = Sphere(Point3f(0), 1.8) |> x -> Tesselation(x, 128) |> uv_mesh

    set_theme!(backgroundcolor=:black)
    scene = Scene(show_axis=false)
    mesh!(
        scene,
        mesh,
        color=texture,
        interpolate=true,
        shading=NoShading
    )
    # point camera to center
    center!(scene)
    cam3d_cad!(scene; fixed_axis=true)
    return scene
end