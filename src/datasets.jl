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
    cells = compute_cell_array(geo_ds.X, geo_ds.Y, resolution, crs)

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

"Fast iterative version only supporting mean"
function to_dggs_dataset(geo_ds::Dataset, resolution::Integer, crs::String; metadata=Dict(), kwargs...)
    cells = compute_cell_array(geo_ds.X, geo_ds.Y, resolution, crs)
    dggs_bbox = get_dggs_bbox(cells)
    geo_bbox = get_geo_bbox(geo_ds.cubes |> values |> first, crs)

    dggs_arrays = []
    Threads.@threads for (name, geo_array) in collect(geo_ds.cubes)
        dggs_array = to_dggs_array(geo_array, cells, dggs_bbox, geo_bbox; name=name, kwargs...)
        push!(dggs_arrays, dggs_array)
    end
    return DGGSDataset(dggs_arrays...; metadata=metadata)
end

function to_geo_dataset(dggs_ds::DGGSDataset, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    cells = compute_cell_array(lon_dim, lat_dim, dggs_ds.resolution)

    geo_arrays = Dict()
    Threads.@threads for k in keys(dggs_ds)
        dggs_array = getproperty(dggs_ds, k)
        geo_array = to_geo_array(dggs_array, cells; kwargs...)
        geo_arrays[k] = geo_array
    end
    geo_ds = Dataset(; geo_arrays...)
    return geo_ds
end

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

open_dggs_dataset(file_path::String; kwargs...) = file_path |> x -> open_dataset(x; kwargs...) |> DGGSDataset

function save_dggs_dataset(file_path::String, dggs_ds::DGGSDataset; kwargs...)
    dggs_ds |> Dataset |> x -> savedataset(x; path=file_path, kwargs...)
end