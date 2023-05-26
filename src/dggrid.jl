using DGGRID7_jll
using GeoJSON
using DataFrames

function dg_call(meta::Dict)
    meta_string = ""
    for (key, val) in meta
        meta_string *= "$(key) $(val)\n"
    end

    tmp_dir = tempname()
    mkdir(tmp_dir)
    meta_path = tempname() # not inside tmp_dir to avoid name collision
    write(meta_path, meta_string)

    DGGRID7_jll.dggrid() do dggrid_path
        cd(tmp_dir)
        run(`$dggrid_path $(meta_path)`)
    end

    rm(meta_path)
    return (tmp_dir)
end

function generate_cells(grid::Grid)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "dggs_type" => "ISEA4H",
        "clip_subset_type" => "WHOLE_EARTH",
        "cell_output_type" => "GEOJSON",
        "cell_output_file_name" => "out"
    )

    if Symbol(grid.type) in GridPresets
        meta["dggs_type"] = grid.type
    else
        meta["dggs_type"] = "CUSTOM"
        meta["dggs_topology"] = grid.topology
        meta["dggs_proj"] = grid.projection
        meta["dggs_res_spec"] = grid.resolution
    end

    dir = dg_call(meta)

    jsonbytes = read("$(dir)/out.geojson")
    fc = GeoJSON.read(jsonbytes, ndim=3)
    df = DataFrame(fc)
    return (df)

    rm(dir, recursive=true)
end