

function Base.show(io::IO, ::MIME"text/plain", cube::CellCube)
    println(io, "DGGS CellCube")
    Base.show(io, "text/plain", cube.data.axes)
end

function CellCube(path::String)
    data = Cube(path)
    level = 6
    CellCube(data, level)
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

function Base.show(io::IO, ::MIME"text/plain", cube::GeoCube)
    println(io, "DGGS GeoCube")
    Base.show(io, "text/plain", cube.data.axes)
end

function Base.getindex(cell_cube::CellCube, i::Q2DI)
    cell_cube.data[q2di_n=At(i.n), q2di_i=At(i.i), q2di_j=At(i.j)]
end

function Base.getindex(cell_cube::CellCube, lon::Real, lat::Real)
    cell_id = transform_points(lon, lat, cell_cube.level)[1, 1]
    cell_cube.data[q2di_n=cell_id.n, q2di_i=At(cell_id.i), q2di_j=At(cell_id.j)]
end

function Base.getindex(cell_cube::CellCube, selector...)
    cell_array = view(cell_cube.data, selector...)
    CellCube(cell_array, cell_cube.data)
end


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

function transform_points(lon_range, lat_range, level)
    points_path = tempname()
    points_string = ""
    # arrange points to match with pixels in png image
    for lon in lon_range
        for lat in lat_range
            points_string *= "$(lon),$(lat)\n"
        end
    end
    write(points_path, points_string)

    out_points_path = tempname()

    meta = Dict(
        "dggrid_operation" => "TRANSFORM_POINTS",
        "dggs_type" => "ISEA4H",
        "dggs_res_spec" => level - 1,
        "input_file_name" => points_path,
        "input_address_type" => "GEO",
        "input_delimiter" => "\",\"", "output_file_name" => out_points_path,
        "output_address_type" => "Q2DI",
        "output_delimiter" => "\",\"",
    )

    call_dggrid(meta)
    cell_ids = CSV.read(out_points_path, DataFrame; header=["q2di_n", "q2di_i", "q2di_j"])
    rm(points_path)
    rm(out_points_path)
    cell_ids_q2di = map((n, i, j) -> Q2DI(n, i, j), cell_ids.q2di_n, cell_ids.q2di_i, cell_ids.q2di_j) |>
                    x -> reshape(x, length(lat_range), length(lon_range))
    return cell_ids_q2di
end

function transform_points(x, y, z, level; tile_length=256)
    # @see https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
    # @see https://help.openstreetmap.org/questions/747/given-a-latlon-how-do-i-find-the-precise-position-on-the-tilew

    longitudes = tile2lng.(range(x, x + 1; length=tile_length), z)
    latitudes = tile2lat.(range(y, y + 1; length=tile_length), z)
    cell_ids = transform_points(longitudes, latitudes, level)
    return cell_ids
end

"Apply function f after filtering of missing and NAN values"
function filter_null(f)
    x -> x |> filter(!ismissing) |> filter(!isnan) |> f
end

function map_geo_to_cell_cube(xout, xin, cell_ids_unique, cell_ids_indexlist, agg_func)
    for (cell_id, cell_indices) in zip(cell_ids_unique, cell_ids_indexlist)
        xout[cell_id.n+1, cell_id.i+1, cell_id.j+1] = agg_func(view(xin, cell_indices))
    end
end

function CellCube(path::String, level; kwargs...)
    geo_cube = GeoCube(path)
    cell_cube = CellCube(geo_cube, 6; kwargs...)
end

function CellCube(geo_cube::GeoCube, level=6, agg_func=filter_null(mean); chunk_size=missing)
    Threads.nthreads() == 1 && @warn "Multithreading is not active. Please consider to start julia with --threads auto"

    @info "Step 1/3: Transform coordinates"

    # precompute spatial mapping (can be reused e.g. for each time point)
    ismissing(chunk_size) ? chunk_size = max(2, 1024 / length(geo_cube.data.lat)) |> ceil |> Int : true

    lon_chunks = Iterators.partition(geo_cube.data.lon, chunk_size) |> collect
    p = Progress(length(lon_chunks))
    cell_ids_mats = @threaded map(lon_chunks) do lons
        next!(p)
        transform_points(lons, geo_cube.data.lat, level)
    end
    finish!(p)

    cell_ids_mat = hcat(cell_ids_mats...) |> permutedims
    cell_ids_unique = unique(cell_ids_mat)

    @info "Step 2/3: Create cell id masks"
    p = Progress(length(cell_ids_unique))
    cell_ids_indexlist = @threaded map(cell_ids_unique) do x
        next!(p)
        findall(isequal(x), cell_ids_mat)
    end
    finish!(p)

    @info "Step 3/3: Re-grid the data"
    cell_cube = mapCube(
        map_geo_to_cell_cube,
        geo_cube.data,
        cell_ids_unique,
        cell_ids_indexlist,
        agg_func,
        indims=InDims(:lon, :lat),
        outdims=OutDims(
            Dim{:q2di_n}(0:11),
            Dim{:q2di_i}(range(0; step=16, length=32)),
            Dim{:q2di_j}(range(0; step=16, length=32)),
        ),
        showprog=true
    )
    return CellCube(cell_cube, level)
end

function map_cell_to_geo_cube(xout, xin, cell_ids_mat, longitudes, latitudes)
    for lon_i in 1:length(longitudes)
        for lat_i in 1:length(latitudes)
            cell_id = cell_ids_mat[lon_i, lat_i]
            xout[lon_i, lat_i] = xin[cell_id.n+1, cell_id.i+1, cell_id.j+1]
        end
    end
end

function GeoCube(cell_cube::CellCube; longitudes=-180:180, latitudes=-90:90)
    # precompute spatial mapping (can be reused e.g. for each time point)
    cell_ids_mat = transform_points(longitudes, latitudes, cell_cube.level)

    geo_array = mapCube(
        map_cell_to_geo_cube,
        cell_cube.data,
        cell_ids_mat,
        longitudes,
        latitudes,
        indims=InDims(:q2di_n, :q2di_i, :q2di_j),
        outdims=OutDims(
            Dim{:lon}(longitudes),
            Dim{:lat}(latitudes)
        )
    )
    return GeoCube(geo_array)
end

function GeoCube(cell_cube::CellCube, x, y, z; cache=missing, tile_length=256)
    # precompute spatial mapping (can be reused e.g. for each time point)
    cell_ids_mat = ismissing(cache) ? transform_points(x, y, z, 6) : cache[x, y, z]
    bbox = BBox(x, y, z)
    longitudes = range(bbox.lat_min, bbox.lat_max; length=tile_length)
    latitudes = range(bbox.lat_min, bbox.lat_max; length=tile_length)
    geo_array = mapCube(
        map_cell_to_geo_cube,
        cell_cube.data,
        cell_ids_mat,
        longitudes,
        latitudes,
        indims=InDims(:q2di_n, :q2di_i, :q2di_j),
        outdims=OutDims(
            Dim{:lon}(longitudes),
            Dim{:lat}(latitudes)
        )
    )
    return GeoCube(geo_array)
end

