using DGGRID7_jll
using NearestNeighbors
using CSV
using DataFrames
using ArchGDAL
using GeoDataFrames
using GeoFormatTypes

struct DgGrid <: AbstractGrid
    data::KDTree
    type::Symbol
    projection::Symbol
    aperture::Int
    topology::Symbol
    level::Int
end

Projections = [:isea, :fuller]
Topologies = [:hexagon, :triangle, :diamond]
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
function get_dggrid_grid_table(type::Symbol, topology::Symbol, projection::Symbol, resolution::Int)
    meta = Dict(
        "dggrid_operation" => "GENERATE_GRID",
        "clip_subset_type" => "WHOLE_EARTH",
        "point_output_type" => "TEXT",
        "point_output_file_name" => "centers"
    )

    meta["dggs_type"] = "CUSTOM"
    meta["dggs_topology"] = uppercase(string(topology))
    meta["dggs_proj"] = uppercase(string(projection))
    meta["dggs_res_spec"] = uppercase(string(resolution))

    out_dir = call_dggrid(meta)

    df = CSV.read("$(out_dir)/centers.txt", DataFrame; header=["name", "lon", "lat"], footerskip=1)
    rm(out_dir, recursive=true)
    return df
end

"""
Create a grid using DGGRID parameters
"""
function DgGrid(projection::Symbol, aperture::Int, topology::Symbol, resolution::Int)
    projection in Projections ? true : error("projection :$projection must be one of $Projections")
    aperture in Apertures ? true : error("aperture $aperture must be one of $Apertures")
    topology in Topologies ? true : error("topology :$(topology) must be one of $Topologies")

    grid_table = get_dggrid_grid_table(:custom, topology, projection, resolution)

    # cell center points encode grid tpopology (e.g. hexagon or square) implicitly
    # Fast average search in O(log n) and efficient in batch processing
    # KDTree defaults to Euklidean metric
    # However, should be faster than haversine and return same indices
    grid_tree = grid_table[:, [:lon, :lat]] |> Matrix |> transpose |> KDTree
    return DgGrid(grid_tree, :custom, projection, aperture, topology, resolution)
end

function Base.show(io::IO, ::MIME"text/plain", grid::DgGrid)
    println(io, "DgGrid with $(grid.topology) topology, $(grid.projection) projection, apterture of $(grid.aperture), and $(length(grid)) cells")
end