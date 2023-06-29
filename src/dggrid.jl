using DGGRID7_jll
import CSV: read
using GeoDataFrames
using DataFrames

"DGGRID projections"
Projections = [:isea, :fuller]

"DGGRID topologies"
Topologies = [:hexagon, :triangle, :diamond]

"DGGRID apertures"
Apertures = [3, 4, 7]

struct Preset
    topology::Symbol
    projection::Symbol
    aperture::Union{Int,Nothing}
end

Presets = Dict(
    :superfund => Preset(:hexagon, :fuller, nothing),
    :planetrisk => Preset(:hexagon, :isea, nothing),
    :isea4t => Preset(:triangle, :isea, 4),
    :isea4d => Preset(:diamond, :isea, 4),
    :isea3h => Preset(:hexagon, :isea, 3),
    :isea4h => Preset(:hexagon, :isea, 4),
    :isea7h => Preset(:hexagon, :isea, 7),
    :isea43h => Preset(:hexagon, :isea, nothing),
    :fuller4t => Preset(:triangle, :fuller, 4),
    :fuller4d => Preset(:diamond, :fuller, 4),
    :fuller3h => Preset(:hexagon, :fuller, 3),
    :fuller4h => Preset(:hexagon, :fuller, 4),
    :fuller7h => Preset(:hexagon, :fuller, 7),
    :fuller43h => Preset(:hexagon, :fuller, nothing)
)

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

"""
Get a DataFrame of cell center points
"""
function get_dggrid_grid_table(topology::Symbol, projection::Symbol, level::Int)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        "point_output_type" => "TEXT",
        "point_output_file_name" => "centers"
    )

    meta["dggs_type"] = "CUSTOM"
    meta["dggs_topology"] = uppercase(string(topology))
    meta["dggs_proj"] = uppercase(string(projection))
    meta["dggs_res_spec"] = uppercase(string(level))

    out_dir = call_dggrid(meta)
    df = read("$(out_dir)/centers.txt", DataFrame; header=["name", "lon", "lat"], footerskip=1)
    rm(out_dir, recursive=true)
    return df
end

function get_dggrid_grid_table(preset::Symbol, level::Int)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        "point_output_type" => "TEXT",
        "point_output_file_name" => "centers"
    )

    meta["dggs_type"] = uppercase(string(preset))
    meta["dggs_res_spec"] = uppercase(string(level))

    out_dir = call_dggrid(meta)
    df = read("$(out_dir)/centers.txt", DataFrame; header=["name", "lon", "lat"], footerskip=1)
    rm(out_dir, recursive=true)
    return df
end

"""
Get a GeoDataFrame with boundary polygons for each cell
"""
function get_dggrid_cell_boundaries(topology::Symbol, projection::Symbol, level::Int)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        # use GDAL geojson, see https://github.com/sahrk/DGGRID/issues/4
        "cell_output_type" => "GDAL",
        "cell_output_gdal_format" => "GeoJSON",
        "cell_output_file_name" => "boundaries.geojson"
    )

    meta["dggs_type"] = "CUSTOM"
    meta["dggs_topology"] = uppercase(string(topology))
    meta["dggs_proj"] = uppercase(string(projection))
    meta["dggs_res_spec"] = uppercase(string(level))

    out_dir = call_dggrid(meta)
    df = GeoDataFrames.read("$(out_dir)/boundaries.geojson")
    rename!(df, [:geometry, :cell_id])
    rm(out_dir, recursive=true)
    return df
end