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
        dggs_array = to_dggs_array(geo_array, cells; kwargs...)
        push!(dggs_arrays, dggs_array)
    end
    return DGGSDataset(dggs_arrays...)
end

function to_geo_dataset(dggs_ds::DGGSDataset, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    cells = compute_cell_array(lon_dim, lat_dim, dggs_ds.resolution)

    geo_arrays = Dict()
    for k in keys(dggs_ds)
        dggs_array = getproperty(dggs_ds, k)
        geo_array = to_geo_array(dggs_array, cells; kwargs...)
        geo_arrays[k] = geo_array
    end
    geo_ds = Dataset(; geo_arrays...)
    return geo_ds
end

Base.show(io::IO, a::DGGSArray) = print(io, "DGGSArray $(name(a)) $(a.dggsrs)@$(a.resolution) ")
Base.show(io::IO, a::DGGSDataset) = print(io, "DGGSDataset $(name(a)) $(a.dggsrs)@$(a.resolution) ")

function print_dggs_block(io, x::Union{DGGSArray,DGGSDataset}; block_width=get(io, :blockwidth, 0))
    DD.print_block_separator(io, "DGGS", block_width, block_width)
    println(io, " ")
    println(io, "  DGGSRS:     $(x.dggsrs)")
    n_cells_str = @sprintf("%.2e", 2 * 2^x.resolution * 2^x.resolution * 5)
    println(io, "  Resolution: $(x.resolution) (up to $(n_cells_str) cells)")
    println(io, "  Geo BBox:   $(get_geo_bbox(x))")
    DD.print_block_close(io, block_width)
end

function DD.show_after(io::IO, mime, x::Union{DGGSArray,DGGSDataset})
    print_dggs_block(io, x)
end