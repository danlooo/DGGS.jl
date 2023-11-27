using Infiltrator

using DGGRID7_jll
using YAXArrays
using DimensionalData
using CSV
using DataFrames
using Statistics
using ColorSchemes

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
        old_pwd = pwd()
        cd(tmp_dir)
        oldstd = stdout
        if !verbose
            redirect_stdout(devnull)
        end
        run(`$dggrid_path $(meta_path)`)
        cd(old_pwd)
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