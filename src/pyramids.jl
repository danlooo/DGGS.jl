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

get_resolutions(dggs_p::DGGSPyramid) = (keys(dggs_p.branches) .|> x -> String(x)[7:end] .|> x -> parse(Int, x)) |> sort

function DD.label(p::DGGSPyramid)
    layer_keys = p.branches |> values |> first |> DD.layers |> keys
    if length(layer_keys) == 1
        return String(layer_keys[1])
    else
        return ""
    end
end

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

function aggregate_by_factor(
    xin::AbstractArray,
    xout::AbstractArray,
    pyramid_agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean
)
    fac = ceil(Int, size(xin, 1) / size(xout, 1))
    for j in axes(xout, 2)
        for i in axes(xout, 1)
            xview = ((i-1)*fac+1):min(size(xin, 1), (i * fac))
            yview = ((j-1)*fac+1):min(size(xin, 2), (j * fac))
            xout[i, j] = pyramid_agg_func(view(xin, xview, yview))
        end
    end
end


function coarsen(
    dggs_array::DGGSArray;
    pyramid_agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean
)
    coarser_level = dggs_array.resolution - 1

    # analog to GeoTIFF overviews 
    coarser_dims = map((:dggs_i, :dggs_j)) do dim
        dim_min, dim_max = dims(dggs_array, dim) |> extrema
        dim_min = floor(dim_min / 2) |> Int
        dim_max = floor(dim_max / 2) |> Int
        Dim{dim}(dim_min:dim_max)
    end

    coarser_arr = mapCube(
        dggs_array;
        indims=InDims(:dggs_i, :dggs_j),
        outdims=OutDims(coarser_dims...)
    ) do xout, xin
        xout = aggregate_by_factor(xin, xout, pyramid_agg_func)
    end

    properties = Dict{String,Any}(metadata(dggs_array))
    properties["dggs_dggsrs"] = dggs_array.dggsrs
    properties["dggs_resolution"] = coarser_level
    properties["dggs_bbox"] = dggs_array.bbox

    coarser_dggs_arr = YAXArray(dims(coarser_arr), coarser_arr.data, properties) |> DGGSArray
    coarser_dggs_arr = rebuild(coarser_dggs_arr; name=name(dggs_array))

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
    data = (pyramid |> reverse .|> x -> x.resolution => x) |> OrderedDict
    pyramid = DGGSPyramid(data, dggs_ds.dggsrs, dggs_ds.bbox)
    return pyramid
end

function to_dggs_pyramid(dggs_array::DGGSArray; kwargs...)
    dggs_ds = DGGSDataset(dggs_array)
    pyramid = to_dggs_pyramid(dggs_ds; kwargs...)
    return pyramid
end

function to_dggs_pyramid(
    geo_ds::YAXArrays.Dataset, resolution::Integer, crs::String, agg_func::Function;
    kwargs...
)
    dggs_ds = to_dggs_dataset(geo_ds, resolution, crs, agg_func; kwargs...)
    dggs_pyramid = to_dggs_pyramid(dggs_ds)
    return dggs_pyramid
end

function to_dggs_pyramid(
    geo_ds::YAXArrays.Dataset,
    resolution::Integer,
    crs::String;
    pyramid_agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean,
    kwargs...
)
    dggs_ds = to_dggs_dataset(geo_ds, resolution, crs; kwargs...)
    dggs_pyramid = to_dggs_pyramid(dggs_ds; pyramid_agg_func=pyramid_agg_func)
    return dggs_pyramid
end

function to_dggs_pyramid(
    geo_array::YAXArrays.YAXArray,
    resolution::Integer,
    crs::String;
    pyramid_agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean,
    kwargs...
)
    dggs_array = to_dggs_array(geo_array, resolution, crs; kwargs...)
    dggs_pyramid = to_dggs_pyramid(dggs_array; pyramid_agg_func=pyramid_agg_func)
    return dggs_pyramid
end

open_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")
save_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")