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
    # Prevent spawning too many threads if called in parallel
    geo_coords = CSV.read(out_points_path, DataFrame; header=["lon", "lat"], ntasks=1)
    geo_coords = map((lon, lat) -> (lon, lat), geo_coords.lon, geo_coords.lat)
    rm(points_path)
    rm(out_points_path)
    return geo_coords
end

"Transforms Vector of (lon,lat) coords to DGGRID indices"
function _transform_points(coords::AbstractVector{Tuple{U,V}}, level) where {U<:Real,V<:Real}
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
    # Prevent spawning too many threads if called in parallel
    cell_ids = CSV.read(out_points_path, DataFrame; header=["q2di_n", "q2di_i", "q2di_j"], ntasks=1)
    rm(points_path)
    rm(out_points_path)
    cell_ids_q2di = map((n, i, j) -> Q2DI(n, i, j), cell_ids.q2di_n, cell_ids.q2di_i, cell_ids.q2di_j)
    return cell_ids_q2di
end

function _transform_points(lon_range, lat_range, level)
    Iterators.product(lon_range, lat_range) |> collect |> vec |> sort |> x -> _transform_points(x, level) |> x -> reshape(x, length(lat_range), length(lon_range))
end

function transform_points(coords::Vector{Tuple{U,V}}, level; show_progress=true, chunk_size_points=2048) where {U<:Real,V<:Real}
    chunks = Iterators.partition(coords, chunk_size_points) |> collect

    if length(chunks) == 1
        return _transform_points(coords, level)
    end

    results = nothing
    if show_progress
        p = Progress(length(chunks))
        results = ThreadsX.map(chunks) do chunk
            next!(p)
            _transform_points(chunk, lat_range, level)
        end
    else
        results = ThreadsX.map(chunks) do chunk
            _transform_points(chunk, lat_range, level)
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
        results = ThreadsX.map(chunks) do chunk
            next!(p)
            _transform_points(chunk, lat_range, level)
        end
    else
        results = ThreadsX.map(chunks) do chunk
            _transform_points(chunk, lat_range, level)
        end
    end

    result = vcat(results...)
    return result
end

"""
chunk_size_points: number of points (e.g. pixels) to transform in one block (task of a thread)
"""
function transform_points(
    lon_range::AbstractVector{A},
    lat_range::AbstractVector{B},
    level::Integer;
    show_progress=true,
    chunk_size_points=2048,
    base_path::String=get(ENV, "DGGS_CELL_IDS_PATH", "https://s3.bgc-jena.mpg.de:9000/dggs/cell_ids")
) where {A<:Real,B<:Real}
    chunk_size_lon = chunk_size_points / length(lat_range) |> ceil |> Int
    lon_chunks = Iterators.partition(lon_range, chunk_size_lon) |> collect

    # try cache
    try
        return open_cell_ids(lon_range, lat_range, level, base_path)
    catch
    end

    # single thread is sufficient
    if length(lon_chunks) == 1
        cell_ids_mat = _transform_points(lon_range, lat_range, level)
        return cell_ids_mat
    end

    cell_ids_mats = nothing
    if show_progress
        p = Progress(length(lon_chunks))
        cell_ids_mats = ThreadsX.map(lon_chunks) do lon_chunk
            next!(p)
            _transform_points(lon_chunk, lat_range, level)
        end
    else
        cell_ids_mats = ThreadsX.map(lon_chunks) do lon_chunk
            _transform_points(lon_chunk, lat_range, level)
        end
    end

    cell_ids_mat = hcat(cell_ids_mats...) |> permutedims
    cell_ids = DimArray(cell_ids_mat, (lon_range |> X, lat_range |> Y))
    return cell_ids
end

function create_cell_ids(
    lon_range::AbstractRange{T},
    lat_range::AbstractRange{T},
    levels::AbstractRange,
    base_path::String=get(ENV, "DGGS_CELL_IDS_PATH", "https://s3.bgc-jena.mpg.de:9000/dggs/cell_ids")
) where {T<:Real}
    for level in levels
        path = "$(base_path)/$(lon_range)_$(lat_range)/$(level)"
        axs = (
            Dim{:lon}(lon_range),
            Dim{:lat}(lat_range),
            Dim{:q2di}(["n", "i", "j"]),
        )
        # avoid complex eltypes e.g (n,i,j) tuple for compatibility reasons and to prevent parsing errors
        a = YAXArray(axs, Zeros(UInt32, length(lon_range), length(lat_range), 3), Dict())
        a = setchunks(a, (2048, 2048, 3)) # do not chunk index dimension
        ds = Dataset(; Dict(:layer => a)...)
        ds = savedataset(ds; path=path, driver=:zarr, skeleton=true, overwrite=true)

        # can not use mapCube here
        # - would require same elementtype of input and output cube
        # - would ignore chunks in slices resulting to cache overflow
        for chunk in ds.layer.chunks
            lon_chunk = chunk[1]
            lat_chunk = chunk[2]

            cell_ids = DGGS.transform_points(lon_range[lon_chunk], lat_range[lat_chunk], level)
            ds.layer[lon=lon_chunk, lat=lat_chunk][:, :, 1] = map(x -> x.n, cell_ids)
            ds.layer[lon=lon_chunk, lat=lat_chunk][:, :, 2] = map(x -> x.i, cell_ids)
            ds.layer[lon=lon_chunk, lat=lat_chunk][:, :, 3] = map(x -> x.j, cell_ids)
        end

        Zarr.consolidate_metadata(path)
    end
end

function open_cell_ids(
    lon_range::AbstractRange{T},
    lat_range::AbstractRange{T},
    level::Integer,
    base_path::String=get(ENV, "DGGS_CELL_IDS_PATH", "https://s3.bgc-jena.mpg.de:9000/dggs/cell_ids")
) where {T<:Real}
    arr = open_dataset("$(base_path)/$(lon_range)_$(lat_range)/$(level)").layer
    # can not set eltype of outdims
    cell_ids = DimArray{Q2DI}(undef, arr.lon, arr.lat)
    mapCube(
        arr,
        indims=InDims(:q2di),
        outdims=OutDims()
    ) do xout, xin
        cell_ids[xin.indices[1], xin.indices[2]] = Q2DI(UInt8(xin[1]), UInt32(xin[2]), UInt32(xin[3]))
    end
    return cell_ids
end