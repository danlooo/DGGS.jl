function DGGSPyramid(data::AbstractDict{T,A}, dggsrs, bbox) where {T,A<:DGGSDataset}
    dimtree = DimTree()
    # add all res levels as branches
    for (resolution, dggs_ds) in pairs(data)
        Base.setproperty!(dimtree, Symbol("dggs_s$(resolution)"), DimTree(dggs_ds))
    end
    res = DGGSPyramid(dimtree, dggsrs, bbox)
    return res
end

function DGGSPyramid(dimtree::DimTree, dggsrs, bbox)
    res = DGGSPyramid(
        getfield(dimtree, :data), dims(dimtree), DD.refdims(dimtree),
        DD.layerdims(dimtree), DD.layermetadata(dimtree), metadata(dimtree),
        DD.branches(dimtree), getfield(dimtree, :tree),
        dggsrs, bbox
    )
    return res
end

Base.propertynames(dggs_p::DGGSPyramid) = union((:dggsrs, :bbox), keys(dggs_p.branches))

function extract_dggs_dataset(dggs_p::DGGSPyramid, layer_name::Symbol)
    # DimTree stores leaves as DimTree objects. Re-create DGGSDataset from layers
    branch = DD.branches(dggs_p)[layer_name]
    dggs_layers = keys(branch)
    arrays = map(x -> branch[x].data, dggs_layers)
    dggs_ds = DGGSDataset(arrays...)
    return dggs_ds
end

function Base.getproperty(p::DGGSPyramid, s::Symbol)
    s in fieldnames(DGGSPyramid) && return getfield(p, s)
    s in keys(p.branches) && return extract_dggs_dataset(p, s)
    error("Key $(s) not found.")
end

function Base.getindex(dggs_p::DGGSPyramid, resolution::Int)
    return getproperty(dggs_p, Symbol("dggs_s$(resolution)"))
end

function DD.show_after(io::IO, mime, x::DGGSPyramid)
    block_width = get(io, :blockwidth, 0)
    DD.print_block_separator(io, "DGGS", block_width, block_width)
    println(io, " ")
    println(io, "  DGGSRS:     $(x.dggsrs)")
    println(io, "  Geo BBox:   $(x.bbox)")
    DD.print_block_close(io, block_width)
end

function aggregate_by_factor(xin::AbstractArray, xout::AbstractArray, f::Function)
    fac = ceil(Int, size(xin, 1) / size(xout, 1))
    for j in axes(xout, 2)
        for i in axes(xout, 1)
            xview = ((i-1)*fac+1):min(size(xin, 1), (i * fac))
            yview = ((j-1)*fac+1):min(size(xin, 2), (j * fac))
            xout[i, j] = f(view(xin, xview, yview))
        end
    end
end

function coarsen(dggs_array::DGGSArray; f=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean)
    coarser_level = dggs_array.resolution - 1

    coarser_arr = mapCube(
        dggs_array;
        indims=InDims(:dggs_i, :dggs_j),
        outdims=OutDims(
            # TODO: restrict range to subset
            Dim{:dggs_i}(range(0; step=1, length=2 * 2^coarser_level)),
            Dim{:dggs_j}(range(0; step=1, length=2^coarser_level))
        )
    ) do xout, xin
        xout = aggregate_by_factor(xin, xout, f)
    end

    properties = Dict{String,Any}(metadata(dggs_array))
    properties["dggs_dggsrs"] = dggs_array.dggsrs
    properties["dggs_resolution"] = coarser_level
    properties["dggs_bbox"] = dggs_array.bbox

    coarser_dggs_arr = YAXArray(dims(coarser_arr), coarser_arr.data, properties) |> DGGSArray
    return coarser_dggs_arr
end

function coarsen(dggs_ds::DGGSDataset; kwargs...)
    coarser_arrays = []
    for key in keys(dggs_ds)
        dggs_array = getproperty(dggs_ds, key)
        coarser_dggs_array = coarsen(dggs_array; kwargs...)
        push!(coarser_arrays, coarser_dggs_array)
    end
    res = DGGSDataset(coarser_arrays...)
    return res
end

function to_dggs_pyramid(dggs_ds::DGGSDataset; kwargs...)
    pyramid = DGGSDataset[]
    push!(pyramid, dggs_ds)
    for resolution in dggs_ds.resolution-1:-1:1
        current_dggs_ds = pyramid[end]
        coarser_ds = coarsen(current_dggs_ds; kwargs...)
        push!(pyramid, coarser_ds)
    end
    data = (pyramid |> reverse .|> x -> x.resolution => x) |> Dict
    pyramid = DGGSPyramid(data, dggs_ds.dggsrs, dggs_ds.bbox)
    return pyramid
end

Base.show(io::IO, p::DGGSPyramid) = print(io, "DGGSPyramid $(p.dggsrs) with resolutions $(first(p.data).second.resolution):$(last(p.data).second.resolution)")

open_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")
save_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")