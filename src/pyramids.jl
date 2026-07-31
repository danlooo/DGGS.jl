"""
    coarsen(A::AbstractArray, factors::Tuple; agg_func=x -> mean(skipmissing(x)))

Coarsen an array by aggregating blocks of elements. Each dimension is reduced
by the corresponding factor.

# Arguments
- `A::AbstractArray`: Input array to coarsen
- `factors::Tuple`: Tuple of coarsening factors, one per dimension.
  Use `1` to keep a dimension unchanged.
- `agg_func`: Aggregation function applied to each block. Default: `mean(skipmissing(x))`.

# Example
```julia
a = rand(64, 32, 10)
coarse_a = coarsen(a, (2, 2, 1))  # Result: 32×16×10
```
"""
function coarsen(A::AbstractArray, new_size::Tuple; agg_func=x -> mean(skipmissing(x)))
    # Build the reshaped dimensions: interleave (new_dim, factor) pairs
    reshaped_dims = Int[]
    for (s, f) in zip(size(A), new_size)
        push!(reshaped_dims, s ÷ f)
        push!(reshaped_dims, f)
    end

    reshaped = reshape(A, Tuple(reshaped_dims))

    # Average over the within-block factor dimensions (every odd dimension: 1, 3, 5, ...)
    # Reshape order is (n_blocks, block_size) per spatial dim. In column-major Julia,
    # n_blocks (even positions) indexes which block, block_size (odd positions) indexes
    # within the block. We reduce over the within-block dimensions.
    reduce_dims = Tuple(1:2:length(reshaped_dims))

    # Compute output shape (keep even dims, drop odd dims)
    out_shape = Tuple(reshaped_dims[i] for i in 2:2:length(reshaped_dims))

    # Determine output eltype by finding first non-missing result
    out_eltype = Missing
    for idx in CartesianIndices(out_shape)
        # Build slice indices: for each reduce_dim, take full range; for output dims, use idx
        slice_indices = []
        out_idx = 1
        for d in 1:length(reshaped_dims)
            if d in reduce_dims
                push!(slice_indices, :)
            else
                push!(slice_indices, idx[out_idx])
                out_idx += 1
            end
        end
        block = reshaped[slice_indices...]
        val = agg_func(block)
        if !ismissing(val)
            out_eltype = typeof(val)
            break
        end
    end

    # Allocate output with correct type
    result = Array{Union{Missing,out_eltype}}(missing, out_shape)

    # Fill output
    for idx in CartesianIndices(out_shape)
        slice_indices = []
        out_idx = 1
        for d in 1:length(reshaped_dims)
            if d in reduce_dims
                push!(slice_indices, :)
            else
                push!(slice_indices, idx[out_idx])
                out_idx += 1
            end
        end
        block = reshaped[slice_indices...]
        result[idx] = agg_func(block)
    end

    return result
end

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

# is defined on data filed by default which is empty in DGGSPyramid
Base.keys(p::DGGSPyramid) = Base.keys(p.branches)
Base.first(p::DGGSPyramid) = p |> keys |> first |> x -> getproperty(p, x)
Base.last(p::DGGSPyramid) = p |> keys |> collect |> last |> x -> getproperty(p, x)

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
    agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean
)
    fac = ceil(Int, size(xin, 1) / size(xout, 1))
    for j in axes(xout, 2)
        for i in axes(xout, 1)
            xview = ((i-1)*fac+1):min(size(xin, 1), (i*fac))
            yview = ((j-1)*fac+1):min(size(xin, 2), (j*fac))
            xout[i, j] = agg_func(view(xin, xview, yview))
        end
    end
end


"""
    coarsen(dggs_array::DGGSArray{<:Any,<:Any,<:Any,<:Any,<:TileArray}; agg_func)

Coarsen a DGGSArray backed by a TileArray by aggregating 2x2 blocks.
Missing tiles are skipped entirely, making this efficient for sparse arrays.
"""
function coarsen(
    dggs_array::DGGSArray{<:Any,<:Any,<:Any,<:Any,<:TileArray};
    agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean
)
    tile_array = dggs_array.data
    coarser_level = dggs_array.resolution - 1

    # Get dimension extents (0-based DGGS coordinates)
    i_min, i_max = extrema(dims(dggs_array, :dggs_i))
    j_min, j_max = extrema(dims(dggs_array, :dggs_j))
    n_min, n_max = extrema(dims(dggs_array, :dggs_n))

    # Compute coarser dimensions (0-based)
    coarser_i_min = floor(Int, i_min / 2)
    coarser_i_max = floor(Int, i_max / 2)
    coarser_j_min = floor(Int, j_min / 2)
    coarser_j_max = floor(Int, j_max / 2)

    # Output TileArray dimensions (1-based extent)
    out_dims = (
        coarser_i_max - coarser_i_min + 1,
        coarser_j_max - coarser_j_min + 1,
        n_max - n_min + 1
    )

    # Create output TileArray with same chunk size
    out_eltype = eltype(tile_array)
    out_tile_array = TileArray{out_eltype}(missing, out_dims, tile_array.chunk_size)

    cs = tile_array.chunk_size

    # Get present tile ranges (1-based internal indices)
    present_ranges = ranges(tile_array)

    # Process each present tile
    for tile_range in present_ranges
        i_range, j_range, n_range = tile_range

        # Get chunk data using 1-based chunk indices
        ci = div(first(i_range) - 1, cs[1]) + 1
        cj = div(first(j_range) - 1, cs[2]) + 1
        cn = div(first(n_range) - 1, cs[3]) + 1
        chunk_data = tile_array.data[ci, cj, cn]
        chunk_data === missing && continue

        # Global 0-based coordinate of chunk start
        global_i_start_0 = first(i_range) - 1 + i_min
        global_j_start_0 = first(j_range) - 1 + j_min

        # Iterate over 2x2 blocks aligned to the global grid
        for li in 1:length(i_range)
            gi_0 = global_i_start_0 + li - 1
            gi_0 % 2 != 0 && continue  # Skip non-aligned positions

            for lj in 1:length(j_range)
                gj_0 = global_j_start_0 + lj - 1
                gj_0 % 2 != 0 && continue

                # Collect values from 2x2 block
                values = out_eltype[]
                for di in 0:1, dj in 0:1
                    ni, nj = li + di, lj + dj
                    if ni <= length(i_range) && nj <= length(j_range)
                        for ln in 1:length(n_range)
                            val = chunk_data[ni, nj, ln]
                            val !== missing && push!(values, val)
                        end
                    end
                end

                if !isempty(values)
                    agg_val = agg_func(values)

                    # Write to output (convert to 1-based)
                    out_i_1 = div(gi_0, 2) - coarser_i_min + 1
                    out_j_1 = div(gj_0, 2) - coarser_j_min + 1
                    for n_idx in n_range
                        out_tile_array[out_i_1, out_j_1, n_idx] = agg_val
                    end
                end
            end
        end
    end

    # Build coarser DGGSArray directly (bypass YAXArray to preserve TileArray)
    coarser_dims = (
        Dim{:dggs_i}(coarser_i_min:coarser_i_max),
        Dim{:dggs_j}(coarser_j_min:coarser_j_max),
        Dim{:dggs_n}(n_min:n_max)
    )

    properties = Dict{String,Any}(metadata(dggs_array))
    properties["dggs_dggsrs"] = dggs_array.dggsrs
    properties["dggs_resolution"] = coarser_level
    properties["dggs_bbox"] = dggs_array.bbox

    coarser_dggs_arr = DGGSArray(
        out_tile_array, coarser_dims, (), name(dggs_array), properties,
        coarser_level, dggs_array.dggsrs, dggs_array.bbox
    )

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
    for resolution in (dggs_ds.resolution-1):-1:1
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
    geo_ds::YAXArrays.Dataset,
    resolution::Integer,
    crs::String;
    agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean,
    kwargs...
)
    dggs_ds = to_dggs_dataset(geo_ds, resolution, crs; agg_func=agg_func, kwargs...)
    dggs_pyramid = to_dggs_pyramid(dggs_ds; agg_func=agg_func)
    return dggs_pyramid
end

function to_dggs_pyramid(
    geo_array::YAXArrays.YAXArray,
    resolution::Integer,
    crs::String;
    agg_func::Function=x -> filter(y -> !ismissing(y) && !isnan(y), x) |> mean,
    kwargs...
)
    dggs_array = to_dggs_array(geo_array, resolution, crs; agg_func=agg_func, kwargs...)
    dggs_pyramid = to_dggs_pyramid(dggs_array; agg_func=agg_func)
    return dggs_pyramid
end

open_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")
save_dggs_pyramid(args...; kwargs...) = error("Please load module Zarr first")