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

function generate_cells(grid::GridSpec)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        "cell_output_type" => "GEOJSON",
        "cell_output_file_name" => "out"
    )

    if Symbol(grid.type) in GridPresets
        meta["dggs_type"] = string(grid.type)
    else
        meta["dggs_type"] = "CUSTOM"
        meta["dggs_topology"] = string(grid.topology)
        meta["dggs_proj"] = string(grid.projection)
        meta["dggs_res_spec"] = string(grid.resolution)
    end

    out_dir = dg_call(meta)

    jsonbytes = read("$(out_dir)/out.geojson")
    cell_feature_collection = GeoJSON.read(jsonbytes, ndim=3)
    df = DataFrame(cell_feature_collection)
    rm(out_dir, recursive=true)
    return (df)
end