using Infiltrator

using DGGRID7_jll
using YAXArrays
using DimensionalData
using CSV
using DataFrames
using Statistics
using ColorSchemes

"Rectangular bounding box in geographical space"
struct BBox{T<:Real}
    lon_min::T
    lon_max::T
    lat_min::T
    lat_max::T
end

"Geographical bounding box given a xyz tile"
function BBox(x, y, z)
    #@see https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames

    n = 2^z
    lon_min = x / n * 360.0 - 180.0
    lat_min = atan(sinh(π * (1 - 2 * (y + 1) / n))) |> rad2deg

    lon_max = (x + 1) / n * 360.0 - 180.0
    lat_max = atan(sinh(π * (1 - 2 * y / n))) |> rad2deg

    return BBox(lon_min, lon_max, lat_min, lat_max)
end

struct ColorScale{T<:Real}
    schema::ColorScheme
    min_value::T
    max_value::T
end

struct Q2DI
    n
    i
    j
end

struct CellCube
    data::YAXArray
    level::Int8
end

struct GeoCube
    data::YAXArray
end

function Base.getindex(cell_cube::CellCube, i::Q2DI)
    cell_cube.data[q2di_n=At(i.n), q2di_i=At(i.i), q2di_j=At(i.j)]
end

function Base.getindex(cell_cube::CellCube, lon::Real, lat::Real)
    cell_id = transform_points(lon, lat, cell_cube.level)[1, 1]
    cell_cube.data[q2di_n=At(cell_id.n), q2di_i=At(cell_id.i), q2di_j=At(cell_id.j)]
end

function Base.getindex(cell_cube::CellCube, selector...)
    cell_array = view(cell_cube.data, selector...)
    CellCube(cell_array, cell_cube.data)
end


"""
Execute sytem call of DGGRID binary
"""
function call_dggrid(meta::Dict; verbose=false)
    meta_string = ""
    for (key, val) in meta
        meta_string *= "$(key) $(val)\n"
    end

    tmp_dir = tempname()
    mkdir(tmp_dir)
    meta_path = tempname() # not inside tmp_dir to avoid name collision
    write(meta_path, meta_string)

    DGGRID7_jll.dggrid() do dggrid_path
        oldstd = stdout
        if !verbose
            redirect_stdout(devnull)
        end

        # ensure thread safetey, e.g. don't use julia functions cd and pwd
        cmd = "cd $tmp_dir && $dggrid_path $meta_path"
        run(`sh -c $cmd`)

        redirect_stdout(oldstd)
    end

    rm(meta_path)
    return (tmp_dir)
end

function transform_points(lon_range, lat_range, level)
    points_path = tempname()
    points_string = ""
    for lat in lat_range
        for lon in lon_range
            points_string *= "$(lon),$(lat)\n"
        end
    end
    write(points_path, points_string)

    meta = Dict(
        "dggrid_operation" => "TRANSFORM_POINTS",
        "dggs_type" => "ISEA4H",
        "dggs_res_spec" => level - 1,
        "input_file_name" => points_path,
        "input_address_type" => "GEO",
        "input_delimiter" => "\",\"", "output_file_name" => "cell_ids.csv",
        "output_address_type" => "Q2DI",
        "output_delimiter" => "\",\"",
    )

    out_dir = call_dggrid(meta)
    cell_ids = CSV.read("$(out_dir)/cell_ids.csv", DataFrame; header=["q2di_n", "q2di_i", "q2di_j"])
    rm(out_dir, recursive=true)
    rm(points_path)
    cell_ids_q2di = map((n, i, j) -> Q2DI(n, i, j), cell_ids.q2di_n, cell_ids.q2di_i, cell_ids.q2di_j) |>
                    x -> reshape(x, length(lon_range), length(lat_range))
    return cell_ids_q2di
end

function transform_points(x, y, z, level; tile_length=256)
    bbox = BBox(x, y, z)
    longitudes = range(bbox.lon_min, bbox.lon_max, tile_length)
    latitudes = range(bbox.lat_min, bbox.lat_max, tile_length)
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

function CellCube(geo_cube::GeoCube, level=6, agg_func=filter_null(mean))
    # precompute spatial mapping (can be reused e.g. for each time point)
    cell_ids_mat = transform_points(geo_cube.data.lon, geo_cube.data.lat, level)
    cell_ids_unique = unique(cell_ids_mat)
    cell_ids_indexlist = map(cell_ids_unique) do x
        findall(isequal(x), cell_ids_mat)
    end

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
        )
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


"""
(pre) calculate x,y,z to cell_ids for lookup cahce in tile server

`max_z`: maximum z level of xyz tiles, result in `max_z` + 1 levels
"""
function calculate_cell_ids_of_tiles(; max_z=3, tile_length=256)
    # flatten tasks to increase multi CPU utilization and 
    tiles_keys = [IterTools.product(0:2^z-1, 0:2^z-1, z) for z in 0:max_z] |> Iterators.flatten |> collect

    # ensure thread saftey. Results might come in differnt order
    result = ThreadSafeDict()
    p = Progress(length(tiles_keys))
    Threads.nthreads() == 1 && @warn "Multithreading is not active. Please consider to start julia with --threads auto"
    Threads.@threads for i in eachindex(tiles_keys)
        tile = tiles_keys[i]
        result[tile...] = DGGS.transform_points(tile[1], tile[2], tile[3], 6)
        next!(p)
    end
    finish!(p)
    return result
end


function color_value(value::Real, color_scale::ColorScale; null_color=RGBA{Float64}(0, 0, 0, 0))
    isnan(value) && return null_color
    ismissing(value) && return null_color
    return color_scale.schema[value] |> RGBA
end


cell_ids = deserialize("cell_ids.jl.dat")

function calculate_tile(x, y, z; tile_length=256)
    color_scale = ColorScale(ColorSchemes.viridis, -4, 4)

    tile_cell_ids = cell_ids[x, y, z]
    tile_values = map(tile_cell_ids) do cell_id
        cell_cube.data[time=2, q2di_n=cell_id.n + 1, q2di_i=cell_id.i + 1, q2di_j=cell_id.j + 1].data[1]
    end

    scaled = (tile_values .- color_scale.min_value) / (color_scale.max_value - color_scale.min_value)
    image = map(x -> color_value(x, color_scale), scaled)

    trfm = LinearMap(RotMatrix2{Float64}(0, -1, 1, 0)) # rotation angle would result in rounding error
    image = warp(image, trfm)

    io = IOBuffer()
    save(Stream(format"PNG", io), image)
    return io.data
end