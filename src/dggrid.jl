
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
