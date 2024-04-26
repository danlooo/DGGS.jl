# Julia - Shell - C++ bridge for DGGRID that implements ISEA grids
# Documentation of DGGRID: https://webpages.sou.edu/~sahrk/docs/dggridManualV70.pdf

"""
Execute sytem call of DGGRID binary
"""
function call_dggrid(meta::Dict)
    meta_string = ""
    for (key, val) in meta
        meta_string *= "$(key) $(val)\n"
    end

    meta_path = tempname()
    write(meta_path, meta_string)

    redirect_stdout(devnull)
    # ensure thread safetey
    # see https://discourse.julialang.org/t/ioerror-could-not-spawn-argument-list-too-long/43728/18
    run(`$(DGGRID7_jll.dggrid()) $meta_path`)

    rm(meta_path)
end

# single threaded version
function _transform_points(coords::AbstractVector{Q2DI{T}}, level) where {T<:Integer}
    points_path = tempname()
    points_string = ""
    for c in coords
        points_string *= "$(c.n), $(c.i), $(c.j)\n"
    end
    write(points_path, points_string)

    out_points_path = tempname()

    meta = Dict(
        "dggrid_operation" => "TRANSFORM_POINTS",
        "dggs_type" => "ISEA4H",
        "dggs_res_spec" => level - 1,
        "input_file_name" => points_path,
        "input_address_type" => "Q2DI",
        "input_delimiter" => "\",\"",
        "output_file_name" => out_points_path,
        "output_address_type" => "GEO",
        "output_delimiter" => "\",\"",
    )

    call_dggrid(meta)
    geo_coords = CSV.read(out_points_path, DataFrame; header=["lon", "lat"])
    geo_coords = map((lon, lat) -> (lon, lat), geo_coords.lon, geo_coords.lat)
    rm(points_path)
    rm(out_points_path)
    return geo_coords
end

"Transforms Vector of (lon,lat) coords to DGGRID indices"
function _transform_points(coords::AbstractVector{Tuple{T,T}}, level) where {T<:Real}
    points_path = tempname()
    points_string = ""
    # arrange points to match with pixels in png image
    for c in coords
        points_string *= "$(c[1]),$(c[2])\n"
    end
    write(points_path, points_string)

    out_points_path = tempname()

    meta = Dict(
        "dggrid_operation" => "TRANSFORM_POINTS",
        "dggs_type" => "ISEA4H",
        "dggs_res_spec" => level - 1,
        "input_file_name" => points_path,
        "input_address_type" => "GEO",
        "input_delimiter" => "\",\"",
        "output_file_name" => out_points_path,
        "output_address_type" => "Q2DI",
        "output_delimiter" => "\",\"",
    )

    call_dggrid(meta)
    cell_ids = CSV.read(out_points_path, DataFrame; header=["q2di_n", "q2di_i", "q2di_j"])
    rm(points_path)
    rm(out_points_path)
    cell_ids_q2di = map((n, i, j) -> Q2DI(n, i, j), cell_ids.q2di_n, cell_ids.q2di_i, cell_ids.q2di_j)
    return cell_ids_q2di
end

function _transform_points(lon_range, lat_range, level)
    product(lon_range, lat_range) |> collect |> vec |> sort |> x -> _transform_points(x, level) |> x -> reshape(x, length(lat_range), length(lon_range))
end

function transform_points(coords::Vector{Tuple{T,T}}, level; show_progress=true, chunk_size_points=2048) where {T<:Real}
    chunks = Iterators.partition(coords, chunk_size_points) |> collect

    if length(chunks) == 1
        return _transform_points(coords, level)
    end

    results = nothing
    if show_progress
        p = Progress(length(chunks))
        results = @threaded map(chunks) do coords
            res = _transform_points(coords, level)
            next!(p)
            res
        end
        finish!(p)
    else
        results = @threaded map(chunks) do coords
            _transform_points(coords, level)
        end
    end

    result = vcat(results...)
    return result
end

function transform_points(coords::Vector{Q2DI{T}}, level; show_progress=true, chunk_size_points=2048) where {T<:Integer}
    chunks = Iterators.partition(coords, chunk_size_points) |> collect

    if length(chunks) == 1
        return _transform_points(coords, level)
    end

    results = nothing
    if show_progress
        p = Progress(length(chunks))
        results = @threaded map(chunks) do coords
            res = _transform_points(coords, level)
            next!(p)
            res
        end
        finish!(p)
    else
        results = @threaded map(chunks) do coords
            _transform_points(coords, level)
        end
    end

    result = vcat(results...)
    return result
end

"""
chunk_size_points: number of points (e.g. pixels) to transform in one block (task of a thread)
"""
function transform_points(lon_range::AbstractVector{A}, lat_range::AbstractVector{B}, level::Integer; show_progress=true, chunk_size_points=2048) where {A<:Real,B<:Real}
    chunk_size_lon = chunk_size_points / length(lat_range) |> ceil |> Int
    lon_chunks = Iterators.partition(lon_range, chunk_size_lon) |> collect

    # single thread is sufficient
    if length(lon_chunks) == 1
        cell_ids_mat = _transform_points(lon_range, lat_range, level)
        return cell_ids_mat
    end

    cell_ids_mats = nothing
    if show_progress
        p = Progress(length(lon_chunks))
        cell_ids_mats = @threaded map(lon_chunks) do lons
            res = _transform_points(lons, lat_range, level)
            next!(p)
            res
        end
        finish!(p)
    else
        cell_ids_mats = @threaded map(lon_chunks) do lons
            _transform_points(lons, lat_range, level)
        end
    end

    cell_ids_mat = hcat(cell_ids_mats...) |> permutedims
    cell_ids = DimArray(cell_ids_mat, (lon_range |> X, lat_range |> Y))
    return cell_ids
end
