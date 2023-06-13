using DGGRID7_jll
using NearestNeighbors
using CSV
using DataFrames
using ArchGDAL
using GeoDataFrames
using GeoFormatTypes

Projections = ["ISEA", "FULLER"]
Topologies = ["HEXAGON", "TRIANGLE", "DIAMOND"]
GridPresets = ["SUPERFUND", "PLANETRISK", "ISEA4T", "ISEA4D", "ISEA3H", "ISEA4H", "ISEA7H", "ISEA43H", "FULLER4T", "FULLER4D", "FULLER3H", "FULLER4H", "FULLER7H", "FULLER43H"]
Apertures = [3, 4, 7]

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
function get_grid_data(grid_spec::GridSpec)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        "point_output_type" => "TEXT",
        "point_output_file_name" => "centers"
    )

    if Symbol(grid_spec.type) in GridPresets
        meta["dggs_type"] = string(grid_spec.type)
    else
        meta["dggs_type"] = "CUSTOM"
        meta["dggs_topology"] = string(grid_spec.topology)
        meta["dggs_proj"] = string(grid_spec.projection)
        meta["dggs_res_spec"] = string(grid_spec.resolution)
    end

    out_dir = call_dggrid(meta)

    df = CSV.read("$(out_dir)/centers.txt", DataFrame; header=["name", "lon", "lat"], footerskip=1)
    rm(out_dir, recursive=true)
    return df
end

get_grid_data(grid::Grid) = get_grid_data(grid.spec)

function get_cell_centers(grid::Grid)
    # Using ArchGDAL directly results in segfaults and code would be more complex
    geometry = Vector{ArchGDAL.IGeometry}(undef, length(grid))
    for i in eachindex(grid.data.data)
        geometry[i] = ArchGDAL.createpoint(grid.data.data[i][1], grid.data.data[i][2])
    end
    df = DataFrame(geometry=geometry, cell_id=grid.data.indices)
    sort!(df, :cell_id)
    return df
end

function export_cell_centers(grid::Grid; filepath::String="centers.geojson")
    df = get_cell_centers(grid)
    GeoDataFrames.write(filepath, df)
end

"""
Get a GeoDataFrame with boundary polygons for each cell
"""
function get_cell_boundaries(grid_spec::GridSpec)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        # use GDAL geojson, see https://github.com/sahrk/DGGRID/issues/4
        "cell_output_type" => "GDAL",
        "cell_output_gdal_format" => "GeoJSON",
        "cell_output_file_name" => "boundaries.geojson"
    )

    if Symbol(grid_spec.type) in GridPresets
        meta["dggs_type"] = string(grid_spec.type)
    else
        meta["dggs_type"] = "CUSTOM"
        meta["dggs_topology"] = string(grid_spec.topology)
        meta["dggs_proj"] = string(grid_spec.projection)
        meta["dggs_res_spec"] = string(grid_spec.resolution)
    end

    out_dir = call_dggrid(meta)
    df = GeoDataFrames.read("$(out_dir)/boundaries.geojson")
    rename!(df, [:geometry, :cell_id])
    rm(out_dir, recursive=true)
    return df
end

get_cell_boundaries(grid::Grid) = get_cell_boundaries(grid.spec)

function export_cell_boundaries(grid::Grid; filepath::String="boundaries.geojson")
    boundaries = get_cell_boundaries(grid)
    GeoDataFrames.write(filepath, boundaries)
end