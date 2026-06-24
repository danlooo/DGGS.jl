"""
ChunkedArray operates just like any AbstractArray, but stores data in chunks. 
This type is optimized for arrays that are dense in only some regions, but most chunks are just not expected to be filled.
In this case, the default value is returned.
This is the in-memory variant of Zarr DictStore, but without Compression and slow hash lookup of the Dict.
Its ideal for global initialized DGGS arrays that only cover a small spatial region, e.g. a couple of UTM tiles.
"""
struct TileArray{T,N} <: AbstractArray{T,N}
    data::Array{Union{Missing,Array{T,N}},N}
    default::T
    dims::NTuple{N,Int}
    chunk_size::NTuple{N,Int}
end

function TileArray{T}(default::T, dims::NTuple{N,Int}, chunk_size::NTuple{N,Int}=dims) where {T,N}
    chunk_dims = ntuple(i -> div(dims[i] + chunk_size[i] - 1, chunk_size[i]), N)
    data = Array{Union{Missing,Array{T,N}},N}(undef, chunk_dims...)
    fill!(data, missing)
    TileArray(data, default, dims, chunk_size)
end

# infer eltype from default if not given
function TileArray(default::T, dims::NTuple{N,Int}, chunk_size::NTuple{N,Int}=dims) where {T,N}
    chunk_dims = ntuple(i -> div(dims[i] + chunk_size[i] - 1, chunk_size[i]), N)
    data = Array{Union{Missing,Array{T,N}},N}(undef, chunk_dims...)
    fill!(data, missing)
    TileArray(data, default, dims, chunk_size)
end

function Base.size(A::TileArray)
    A.dims
end

function Base.length(A::TileArray)
    prod(A.dims)
end

function Base.eltype(::TileArray{T,N}) where {T,N}
    return T
end

Base.IndexStyle(::Type{<:TileArray}) = IndexCartesian()

function Base.getindex(A::TileArray, I::Vararg{Int,N}) where {N}
    if length(I) != length(A.dims)
        throw(DimensionMismatch("Number of indices does not match array dimensions"))
    end
    chunk_key = ntuple(i -> div(I[i] - 1, A.chunk_size[i]) + 1, N)
    chunk = A.data[chunk_key...]
    if chunk === missing
        return A.default
    else
        local_indices = ntuple(i -> mod(I[i] - 1, A.chunk_size[i]) + 1, N)
        return chunk[local_indices...]
    end
end

function Base.setindex!(A::TileArray, value, I::Vararg{Int,N}) where {N}
    chunk_key = ntuple(i -> div(I[i] - 1, A.chunk_size[i]) + 1, N)
    if A.data[chunk_key...] === missing
        A.data[chunk_key...] = fill(A.default, A.chunk_size...)
    end
    local_indices = ntuple(i -> mod(I[i] - 1, A.chunk_size[i]) + 1, N)
    A.data[chunk_key...][local_indices...] = value
end

function Base.iterate(A::TileArray, state=1)
    if state > length(A)
        return nothing
    end
    I = ntuple(i -> div(state - 1, prod(A.dims[(i+1):end])) % A.dims[i] + 1, length(A.dims))
    return (getindex(A, I...), state + 1)
end

function Base.show(io::IO, ::MIME"text/plain", a::TileArray)
    print(io, join(a.dims, "x"))
    print(io, " ")
    print(io, typeof(a))
end

"""
    _tile_reduce_contribution(op, base_result, n::Int)

Given `base_result = f(x)` and the fact that value `x` appears `n` times,
return the reduction of these `n` identical values: `reduce(op, fill(base_result, n))`.
"""
function _tile_reduce_contribution(op, base_result, n::Int)
    n <= 0 && throw(DomainError(n, "n must be positive"))
    n == 1 && return base_result
    if op === +
        return n * base_result
    end
    result = base_result
    remaining = n - 1
    while remaining > 0
        result = op(result, base_result)
        remaining -= 1
    end
    result
end

# Actual size of chunk `ck_d` (1-based) along dimension `d` of array `A`,
# accounting for boundary chunks that may be smaller than `chunk_size[d]`.
@inline _tile_actual_chunk_dim(A, d, ck_d::Int) =
    min(A.chunk_size[d], A.dims[d] - (ck_d - 1) * A.chunk_size[d])

"""
    Base.mapreduce(f, op, A::TileArray; dims=:, init)

Efficiently compute `mapreduce(f, op, A)` over a `TileArray`.

Missing tiles (which contain only `A.default`) are never materialized:
their contribution is computed analytically from the count of missing
elements and the value `f(A.default)`.

Supported forms:
- `mapreduce(f, op, A::TileArray)` — full reduction to a scalar.
- `mapreduce(f, op, A::TileArray; dims=...)` — reduction along `dims`,
  returning an array of shape `(d -> d in dims ? 1 : size(A, d))`.
- `init` optionally seeds the reduction; with `init=Val(:init)` a value can be
  passed as the fourth positional argument following Base's convention.
"""
function Base.mapreduce(f, op, A::TileArray{T,N}; dims=:, init=nothing) where {T,N}
    # When default is missing, missing tiles are "empty" and should be skipped
    A.default === missing && return _mapreduce_tile_skipmissing(A, f, op, dims, init)
    _mapreduce_tile(A, f, op, dims, init)
end

# Internal: mapreduce for TileArray where A.default === missing
# Skip all missing tiles entirely and treat missing elements within chunks as absent
function _mapreduce_tile_skipmissing(A::TileArray{T,N}, f, op, dims, init) where {T,N}
    chunk_dims_tuple = size(A.data)

    if dims === (:)
        # -------------- Full reduction to scalar --------------
        acc = nothing
        for ck in CartesianIndices(chunk_dims_tuple)
            chunk = A.data[ck]
            chunk === missing && continue  # skip missing tiles
            for v in chunk
                v === missing && continue  # skip missing elements within chunks
                val = f(v)
                acc = acc === nothing ? val : op(acc, val)
            end
        end
        if acc === nothing
            return init === nothing ?
                   throw(ArgumentError("reducing over an empty collection with no initial value is not allowed")) :
                   init
        end
        return init === nothing ? acc : op(init, acc)
    end

    # -------------- Reduction along `dims` to an array --------------
    dims_tuple = dims isa Union{Tuple,AbstractVector} ? Tuple(dims) : (dims,)
    out_dims = ntuple(d -> d in dims_tuple ? 1 : A.dims[d], N)
    OutEltype = if init === nothing
        typeof(first(_nonmiss_from_present_chunk(A)))
    else
        typeof(op(init, first(_nonmiss_from_present_chunk(A))))
    end

    present_acc = Array{Union{Nothing,OutEltype},N}(undef, out_dims...)
    fill!(present_acc, nothing)
    seen_present = fill(false, out_dims...)

    for ck in CartesianIndices(chunk_dims_tuple)
        chunk = A.data[ck]
        chunk === missing && continue
        local_shape = ntuple(d -> _tile_actual_chunk_dim(A, d, ck[d]), N)
        for li in CartesianIndices(local_shape)
            v = chunk[li]
            v === missing && continue
            li_vals = Tuple(li)
            gi = ntuple(d -> (ck[d] - 1) * A.chunk_size[d] + li_vals[d], N)
            out_idx = ntuple(d -> d in dims_tuple ? 1 : gi[d], N)
            val = f(v)
            if !seen_present[out_idx...]
                present_acc[out_idx...] = val
                seen_present[out_idx...] = true
            else
                present_acc[out_idx...] = op(present_acc[out_idx...], val)
            end
        end
    end

    acc = Array{OutEltype,N}(undef, out_dims...)
    for out_idx in CartesianIndices(out_dims...)
        p = present_acc[out_idx]
        if init === nothing
            if p === nothing
                acc[out_idx] = zero(OutEltype)
            else
                acc[out_idx] = p
            end
        else
            acc[out_idx] = p !== nothing ? op(init, p) : init
        end
    end
    return acc
end

# Handle `mapreduce` over `skipmissing(a)`. Missing tiles (chunks of only `missing`)
# are skipped entirely — we only fold over elements from present chunks.
function Base.mapreduce(f, op, itr::Base.SkipMissing{<:TileArray{T,N}}; dims=:, init=nothing) where {T,N}
    A = itr.x
    chunk_dims_tuple = size(A.data)

    if dims === (:)
        acc = nothing
        for ck in CartesianIndices(chunk_dims_tuple)
            chunk = A.data[ck]
            chunk === missing && continue  # whole tile is missing — skip
            for v in chunk
                v === missing && continue
                val = f(v)
                acc = acc === nothing ? val : op(acc, val)
            end
        end
        if acc === nothing
            return init === nothing ?
                   throw(ArgumentError("reducing over an empty collection with no initial value is not allowed")) :
                   init
        end
        return init === nothing ? acc : op(init, acc)
    end

    # Reduction along dims — delegate to generic TileArray reduction that skips missing values
    dims_tuple = dims isa Union{Tuple,AbstractVector} ? Tuple(dims) : (dims,)
    out_dims = ntuple(d -> d in dims_tuple ? 1 : A.dims[d], N)
    OutEltype = if init === nothing
        typeof(first(_nonmiss_from_present_chunk(A)))
    else
        typeof(op(init, first(_nonmiss_from_present_chunk(A))))
    end
    miss_count = zeros(Int, out_dims...)  # not used for skipmissing but kept for uniformity
    present_acc = Array{Union{Nothing,OutEltype},N}(undef, out_dims...)
    fill!(present_acc, nothing)
    seen_present = fill(false, out_dims...)

    for ck in CartesianIndices(chunk_dims_tuple)
        chunk = A.data[ck]
        chunk === missing && continue
        local_shape = ntuple(d -> _tile_actual_chunk_dim(A, d, ck[d]), N)
        for li in CartesianIndices(local_shape)
            v = chunk[li]
            v === missing && continue
            li_vals = Tuple(li)
            gi = ntuple(d -> (ck[d] - 1) * A.chunk_size[d] + li_vals[d], N)
            out_idx = ntuple(d -> d in dims_tuple ? 1 : gi[d], N)
            val = f(v)
            if !seen_present[out_idx...]
                present_acc[out_idx...] = val
                seen_present[out_idx...] = true
            else
                present_acc[out_idx...] = op(present_acc[out_idx...], val)
            end
        end
    end

    acc = Array{OutEltype,N}(undef, out_dims...)
    for out_idx in CartesianIndices(out_dims...)
        p = present_acc[out_idx]
        if init === nothing
            if p === nothing
                acc[out_idx] = zero(OutEltype)  # shouldn't happen for skipmissing with data
            else
                acc[out_idx] = p
            end
        else
            acc[out_idx] = p !== nothing ? op(init, p) : init
        end
    end
    return acc
end

# Helper: find first non-missing element from any present chunk (used only for type inference)
function _nonmiss_from_present_chunk(A::TileArray{T,N}) where {T,N}
    for chunk in A.data
        chunk === missing && continue
        for v in chunk
            v === missing || return Some(v)
        end
    end
    # Fallback if all present chunks are missing (shouldn't happen in practice)
    return Some(A.default)
end
Base.first(s::Some) = s.value

# Internal entry point. The `kw...` splat lets us also accept the Base-style
# four-argument form `mapreduce(f, op, ::InitialValueOperator, A)` but here
# we only handle the keyword-based `init`.
function _mapreduce_tile(A::TileArray{T,N}, f, op, dims, init) where {T,N}
    f_default = f(A.default)
    chunk_dims_tuple = size(A.data)
    skip_missing = f_default === missing

    if dims === (:)
        # -------------- Full reduction to scalar --------------
        present_acc = nothing
        n_missing_elems = 0
        for ck in CartesianIndices(chunk_dims_tuple)
            chunk = A.data[ck]
            if chunk === missing
                n_missing_elems += prod(_tile_actual_chunk_dim(A, d, ck[d]) for d in 1:N)
            else
                if skip_missing
                    # Skip missing elements within the chunk
                    for v in chunk
                        v === missing && continue
                        val = f(v)
                        present_acc = present_acc === nothing ? val : op(present_acc, val)
                    end
                else
                    chunk_val = mapreduce(f, op, chunk)
                    present_acc = present_acc === nothing ? chunk_val : op(present_acc, chunk_val)
                end
            end
        end
        result = if n_missing_elems > 0 && !skip_missing
            mc = _tile_reduce_contribution(op, f_default, n_missing_elems)
            present_acc === nothing ? mc : op(present_acc, mc)
        else
            present_acc
        end
        if result === nothing
            if init === nothing
                throw(ArgumentError(
                    "reducing over an empty collection with no initial value is not allowed"))
            end
            return init
        end
        return init === nothing ? result : op(init, result)
    end

    # -------------- Reduction along `dims` to an array --------------
    dims_tuple = dims isa Union{Tuple,AbstractVector} ? Tuple(dims) : (dims,)
    out_dims = ntuple(d -> d in dims_tuple ? 1 : A.dims[d], N)

    # Infer the output element type.
    OutEltype = if init === nothing
        typeof(f_default)
    else
        typeof(op(init, f_default))
    end

    # Accumulators: per output cell, store (i) count of `default` elements,
    # (ii) reduction over values from present chunks. The latter is `nothing`
    # until at least one present element folds into it.
    miss_count = zeros(Int, out_dims...)
    present_acc = Array{Union{Nothing,OutEltype},N}(undef, out_dims...)
    fill!(present_acc, nothing)

    # Scatter pass: visit each chunk once.
    # For a missing chunk covering non-reduced element indices `I_nr`, each
    # output cell `out_idx` indexed by `I_nr` gets `n_miss = prod(actual_size
    # along reduced dims)` elements of `A.default`. For a present chunk, we
    # iterate every element and fold `f(chunk[I])` into the corresponding
    # output cell.
    #
    # For missing chunks that span many non-reduced output cells, this is
    # cheaper than the gather approach because we don't re-iterate the
    # reduced-dim chunk range for each output cell.
    for ck in CartesianIndices(size(A.data))
        chunk = A.data[ck]
        # Reduced-dim actual sizes of this chunk:
        n_miss = prod(_tile_actual_chunk_dim(A, d, ck[d]) for d in dims_tuple)
        if chunk === missing
            # Iterate over non-reduced element indices within this chunk.
            # For each such index, the output-cell index has non-reduced coords
            # equal to the global element index and reduced coords equal to 1.
            non_red_ranges = ntuple(d -> d in dims_tuple ? (1:1) :
                                         (1:_tile_actual_chunk_dim(A, d, ck[d])), N)
            for li_nr in CartesianIndices(non_red_ranges)
                li_vals = Tuple(li_nr)
                out_idx = ntuple(d -> d in dims_tuple ? 1 :
                                      (ck[d] - 1) * A.chunk_size[d] + li_vals[d], N)
                @inbounds miss_count[out_idx...] += n_miss
            end
        else
            # Present chunk: fold contributions into per-output-cell accumulator.
            local_shape = ntuple(d -> _tile_actual_chunk_dim(A, d, ck[d]), N)
            for li in CartesianIndices(local_shape)
                li_vals = Tuple(li)
                gi = ntuple(d -> (ck[d] - 1) * A.chunk_size[d] + li_vals[d], N)
                out_idx = ntuple(d -> d in dims_tuple ? 1 : gi[d], N)
                val = f(chunk[li])
                prev = present_acc[out_idx...]
                @inbounds present_acc[out_idx...] = prev === nothing ? val : op(prev, val)
            end
        end
    end

    # Final pass: merge miss_count and present_acc into the output array.
    # Skip missing-tile contributions when f_default === missing (they're "empty" tiles).
    skip_missing_default = f_default === missing
    acc = Array{OutEltype,N}(undef, out_dims...)
    for out_idx in CartesianIndices(out_dims)
        m = miss_count[out_idx]
        p = present_acc[out_idx]
        if skip_missing_default
            m = 0  # treat missing-tile elements as absent
        end
        if init === nothing
            if m == 0 && p === nothing
                throw(ArgumentError(
                    "reducing over an empty collection with no initial value is not allowed"))
            elseif m == 0
                acc[out_idx] = p
            elseif p === nothing
                acc[out_idx] = _tile_reduce_contribution(op, f_default, m)
            else
                acc[out_idx] = op(_tile_reduce_contribution(op, f_default, m), p)
            end
        else
            base = init
            if m > 0
                base = op(base, _tile_reduce_contribution(op, f_default, m))
            end
            if p !== nothing
                base = op(base, p)
            end
            acc[out_idx] = base
        end
    end
    return acc
end