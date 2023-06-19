using DGGRID7_jll
using NearestNeighbors
using CSV
using DataFrames
using ArchGDAL
using GeoDataFrames
using GeoFormatTypes

struct DgGrid <: AbstractGrid
    data::KDTree
    type::String
    projection::String
    aperture::Int
    topology::String
    resolution::Int
end


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
function get_dggrid_grid_table(type::String, topology::String, projection::String, resolution::Int)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        "point_output_type" => "TEXT",
        "point_output_file_name" => "centers"
    )

    if Symbol(type) in GridPresets
        meta["dggs_type"] = string(type)
    else
        meta["dggs_type"] = "CUSTOM"
        meta["dggs_topology"] = string(topology)
        meta["dggs_proj"] = string(projection)
        meta["dggs_res_spec"] = string(resolution)
    end

    out_dir = call_dggrid(meta)

    df = CSV.read("$(out_dir)/centers.txt", DataFrame; header=["name", "lon", "lat"], footerskip=1)
    rm(out_dir, recursive=true)
    return df
end

"""
Create a grid using DGGRID parameters
"""
function DgGrid(projection::String, aperture::Int, topology::String, resolution::Int)
    if !(projection in Projections)
        throw(DomainError("Argument projection must be an any of $(join(Projections, ","))"))
    end

    if !(aperture in Apertures)
        throw(DomainError("Argument aperture must be an any of $(join(Apertures, ","))"))
    end

    if !(topology in Topologies)
        throw(DomainError("Argument topology must be an any of $(join(Topologies, ","))"))
    end

    grid_table = get_dggrid_grid_table("CUSTOM", topology, projection, resolution)

    # cell center points encode grid tpopology (e.g. hexagon or square) implicitly
    # Fast average search in O(log n) and efficient in batch processing
    # KDTree defaults to Euklidean metric
    # However, should be faster than haversine and return same indices
    grid_tree = grid_table[:, [:lon, :lat]] |> Matrix |> transpose |> KDTree
    return DgGrid(grid_tree, "CUSTOM", projection, aperture, topology, resolution)
end
