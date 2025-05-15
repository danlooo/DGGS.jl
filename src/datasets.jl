function DGGSDataset(dggs_array::DGGSArray)
    ds = DimStack(dggs_array)
    arrays = (; DD.name(dggs_array) => dggs_array)

    return DGGSDataset(
        arrays, dims(ds), DD.refdims(ds), DD.layerdims(ds), metadata(ds),
        DD.layermetadata(ds), dggs_array.resolution, dggs_array.dggsrs
    )
end

function DGGSDataset(dggs_arrays...)
    all(map(x -> x isa DGGSArray, dggs_arrays)) || error("All arrays must be of type DGGSArray")

    resolution = dggs_arrays[1].resolution
    all(map(x -> x.resolution, dggs_arrays) .== resolution) || error("Resolution of all DGGS arrays must be the same")

    dggsrs = dggs_arrays[1].dggsrs
    all(map(x -> x.dggsrs, dggs_arrays) .== dggsrs) || error("DGGSRS of all DGGS arrays must be the same")

    if map(x -> DD.name(x), dggs_arrays) |> unique |> length != length(dggs_arrays)
        error("Name of DGGS arrays must be unique")
    end

    array_tuple = map(x -> DD.name(x) => x, dggs_arrays) |> NamedTuple
    ds = DimStack(array_tuple)

    return DGGSDataset(
        array_tuple, dims(ds), DD.refdims(ds), DD.layerdims(ds), metadata(ds),
        DD.layermetadata(ds), resolution, dggsrs
    )
end

Base.propertynames(ds::DGGSDataset) = union((:resolution, :dggsrs), Base.propertynames(parent(ds)))

function Base.getproperty(ds::DGGSDataset, s::Symbol)
    if s in keys(ds)
        return DGGSArray(ds[s], ds.resolution, ds.dggsrs)
    end
    s in fieldnames(DGGSDataset) && return getfield(ds, s)
    error("Key $(s) not found.")
end

function to_dggs_dataset(geo_ds::Dataset, resolution::Integer, crs::String; kwargs...)
    cells = compute_cell_array(geo_ds.X, geo_ds.Y, resolution, crs)
    dggs_arrays = []
    for (name, geo_array) in geo_ds.cubes
        dggs_array = to_dggs_array(geo_array, cells, geo_ds.X, geo_ds.Y; kwargs...)
        push!(dggs_arrays, dggs_array)
    end
    return DGGSDataset(dggs_arrays...)
end

function DD.show_after(io::IO, mime, ds::DGGSDataset)
    block_width = get(io, :blockwidth, 0)

    DD.print_block_separator(io, "DGGS", block_width, block_width)
    println(io, "  DGGSRS:     $(ds.dggsrs)")
    println(io, "  Resolution: $(ds.resolution)")
    DD.print_block_close(io, block_width)
end

Base.show(io::IO, a::DGGSArray) = print(io, "DGGSArray $(name(a)) $(a.dggsrs)@$(a.resolution) ")
Base.show(io::IO, a::DGGSDataset) = print(io, "DGGSDataset $(name(a)) $(a.dggsrs)@$(a.resolution) ")

function print_dggs_block(io, x::Union{DGGSArray,DGGSDataset}; block_width=get(io, :blockwidth, 0))
    DD.print_block_separator(io, "DGGS", block_width, block_width)
    println(io, " ")
    println(io, "  DGGSRS:     $(x.dggsrs)")
    println(io, "  Resolution: $(x.resolution)")
    DD.print_block_close(io, block_width)
end

function DD.show_after(io::IO, mime, x::Union{DGGSArray,DGGSDataset})
    print_dggs_block(io, x)
end