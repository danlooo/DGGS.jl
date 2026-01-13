function DGGSDataset(dggs_arrays...; kwargs...)
    all(map(x -> x isa DGGSArray, dggs_arrays)) || error("All arrays must be of type DGGSArray")

    resolution = dggs_arrays[1].resolution
    all(map(x -> x.resolution, dggs_arrays) .== resolution) || error("Resolution of all DGGS arrays must be the same")

    dggsrs = dggs_arrays[1].dggsrs
    all(map(x -> x.dggsrs, dggs_arrays) .== dggsrs) || error("DGGSRS of all DGGS arrays must be the same")

    bbox = dggs_arrays[1].bbox
    all(dggs_arrays .|> x -> isequal(x.bbox, bbox)) || error("Bounding box of all DGGS arrays must be the same")

    if map(x -> DD.name(x), dggs_arrays) |> unique |> length != length(dggs_arrays)
        error("Name of DGGS arrays must be unique")
    end

    extract_name(a) = a.name == DD.NoName() ? :layer1 : a.name
    array_tuple = map(x -> extract_name(x) => x, dggs_arrays) |> NamedTuple
    ds = DimStack(array_tuple; kwargs...)

    return DGGSDataset(
        array_tuple, dims(ds), DD.refdims(ds), DD.layerdims(ds), metadata(ds),
        DD.layermetadata(ds), resolution, dggsrs, bbox
    )
end

Base.propertynames(ds::DGGSDataset) = union((:resolution, :dggsrs), Base.propertynames(parent(ds)))

function Base.getproperty(ds::DGGSDataset, s::Symbol)
    s in keys(ds) && return ds.data[s]
    s in fieldnames(DGGSDataset) && return getfield(ds, s)
    error("Key $(s) not found.")
end

function to_dggs_dataset(geo_ds::Dataset, resolution::Integer, crs::String, agg_func::Function; metadata=Dict(), kwargs...)
    cells = to_cell_array(geo_ds.X, geo_ds.Y, resolution, crs)

    # get pixels to aggregate for each cell
    cell_coords = Dict{eltype(cells),Vector{CartesianIndex{2}}}()
    for cell_idx in CartesianIndices(cells)
        cell = cells[cell_idx]
        current_cells = get!(() -> CartesianIndex{2}[], cell_coords, cell)
        push!(current_cells, cell_idx)
    end

    dggs_bbox = get_dggs_bbox(keys(cell_coords))

    dggs_arrays = []
    Threads.@threads for (name, geo_array) in collect(geo_ds.cubes)
        dggs_array = to_dggs_array(geo_array, cells, cell_coords, dggs_bbox, agg_func; name=name, kwargs...)
        push!(dggs_arrays, dggs_array)
    end
    return DGGSDataset(dggs_arrays...; metadata=metadata)
end

"Fast iterative version using nearest neioghbor pixel"
function to_dggs_dataset(
    geo_ds::Dataset,
    resolution::Integer,
    crs::String;
    path=tempname() * ".dggs.zarr",
    x_name=:X,
    y_name=:Y,
    metadata=Dict(),
    n_parallel_chunks=Threads.nthreads(),
    chunks=(dggs_i=4096, dggs_j=4096, dggs_n=1),
    kwargs...
)
    geo_ds = cache(geo_ds)
    dggs_ds = DGGS.init_global_dggs_dataset(geo_ds, resolution, crs, path; kwargs...)

    x_dim = geo_ds.axes[x_name] |> DD.format
    y_dim = geo_ds.axes[y_name] |> DD.format
    used_chunks = get_chunks(resolution, x_dim, y_dim, crs; chunks=chunks)
    batches = Iterators.partition(used_chunks, n_parallel_chunks)

    # limit memory usage by processing in batches
    for batch in batches
        # multi-threading without task migration making proj thread safe
        Threads.@threads :static for chunk in batch
            # pre-allocate chunk data
            chunk_arrays = Dict()
            for (array_name, geo_array) in pairs(geo_ds.cubes)
                non_spatial_dims = filter(x -> !(name(x) in [x_name, y_name]), dims(geo_array))
                chunk_spatial_dims = (
                    dims(dggs_ds, :dggs_i)[chunk[1]],
                    dims(dggs_ds, :dggs_j)[chunk[2]],
                    dims(dggs_ds, :dggs_n)[chunk[3]]
                )
                chunk_dims = (chunk_spatial_dims..., non_spatial_dims...)
                chunk_lengths = length.(chunk_dims)
                chunk_data = Array{eltype(geo_array)}(undef, chunk_lengths...)
                chunk_array = DimArray(chunk_data, chunk_dims)
                chunk_arrays[array_name] = chunk_array
            end

            # cache for other dims, e.g. time steps and bands
            trans = Proj.Transformation(DGGS.crs_isea, crs, ctx=Proj.proj_context_create(), always_xy=true)
            for (i, j, n) in Iterators.product(chunk...)
                cell = Cell(i - 1, j - 1, n - 1, resolution)
                x, y = to_geo(cell, trans)

                point_geo_ds = geo_ds[x=Near(x), Y=Near(y)]
                for (array_name, geo_array) in pairs(point_geo_ds.cubes)
                    chunk_data = length(geo_array) == 1 ? geo_array[1] : geo_array
                    chunk_arrays[array_name][At(cell.i), At(cell.j), At(cell.n), :] = chunk_data
                end
            end

            @sync for (array_name, dggs_array) in pairs(dggs_ds.data)
                Threads.@spawn begin
                    # write in chunks
                    dggs_ds[array_name][dggs_i=chunk[1], dggs_j=chunk[2], dggs_n=chunk[3]] = chunk_arrays[array_name]
                end
            end
        end
    end

    return dggs_ds
end

function to_geo_dataset(dggs_ds::DGGSDataset, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    cells = to_cell_array(lon_dim, lat_dim, dggs_ds.resolution)

    geo_arrays = Dict()
    Threads.@threads for k in keys(dggs_ds)
        dggs_array = getproperty(dggs_ds, k)
        geo_array = to_geo_array(dggs_array, cells; kwargs...)
        geo_arrays[k] = geo_array
    end
    geo_ds = Dataset(; geo_arrays...)
    return geo_ds
end

init_global_dggs_dataset() = @error("Please load package Zarr")

Base.show(io::IO, a::DGGSArray) = print(io, "DGGSArray $(name(a)) $(a.dggsrs)@$(a.resolution) ")
Base.show(io::IO, a::DGGSDataset) = print(io, "DGGSDataset $(name(a)) $(a.dggsrs)@$(a.resolution) ")

function DD.show_after(io::IO, mime, x::Union{DGGSArray,DGGSDataset})
    block_width = get(io, :blockwidth, 0)
    DD.print_block_separator(io, "DGGS", block_width, block_width)
    println(io, " ")
    println(io, "  DGGSRS:     $(x.dggsrs)")
    n_cells_str = @sprintf("%.2e", 2 * 2^x.resolution * 2^x.resolution * 5)
    println(io, "  Resolution: $(x.resolution) (up to $(n_cells_str) cells)")
    println(io, "  Geo BBox:   $(x.bbox)")
    DD.print_block_close(io, block_width)
end

#
# IO:: Serialization of DGGS Datasets
#

DGGSDataset(ds::Dataset) = DGGSDataset([DGGSArray(v) for (k, v) in ds.cubes]...; metadata=ds.properties)

function YAXArrays.Dataset(dggs_ds::DGGSDataset)
    arrays = [k => getproperty(dggs_ds, k) |> YAXArray for k in keys(dggs_ds)]
    properties = Dict{String,Any}(metadata(dggs_ds))
    properties["dggs_dggsrs"] = dggs_ds.dggsrs
    properties["dggs_resolution"] = dggs_ds.resolution
    Dataset(; properties=properties, arrays...)
end

open_dggs_dataset(file_path::String; kwargs...) = file_path |> x -> open_dataset(x; kwargs...) |> cache |> DGGSDataset

function save_dggs_dataset(file_path::String, dggs_ds::DGGSDataset; kwargs...)
    dggs_ds |> Dataset |> x -> savedataset(x; path=file_path, kwargs...)
end