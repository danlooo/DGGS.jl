"""
ChunkedArray operates just like any AbstractArray, but stores data in chunks. 
This type is optimized for arrays that are dense in only some regions, but most chunks are just not expected to be filled.
In this case, the default value is returned.
This is the in-memory variant of Zarr DictStore, but without Compression and slow hash lookup of the Dict.
Its ideal for global initialized DGGS arrays that only cover a small spatial region, e.g. a couple of UTM tiles.
"""
struct ChunkedArray{T,N} <: AbstractArray{T,N}
    data::Array{Union{Missing,Array{T,N}},N}
    default::T
    dims::NTuple{N,Int}
    chunk_size::NTuple{N,Int}
end

function ChunkedArray{T}(default::T, dims::NTuple{N,Int}, chunk_size::NTuple{N,Int}=dims) where {T,N}
    chunk_dims = ntuple(i -> div(dims[i] + chunk_size[i] - 1, chunk_size[i]), N)
    data = Array{Union{Missing,Array{T,N}},N}(undef, chunk_dims...)
    fill!(data, missing)
    ChunkedArray(data, default, dims, chunk_size)
end

function Base.size(A::ChunkedArray)
    A.dims
end

function Base.length(A::ChunkedArray)
    prod(A.dims)
end

function Base.eltype(A::ChunkedArray)
    typeof(A.default)
end

Base.IndexStyle(::Type{<:ChunkedArray}) = IndexCartesian()

function Base.getindex(A::ChunkedArray, I::Vararg{Int,N}) where {N}
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

function Base.setindex!(A::ChunkedArray, value, I::Vararg{Int,N}) where {N}
    chunk_key = ntuple(i -> div(I[i] - 1, A.chunk_size[i]) + 1, N)
    if A.data[chunk_key...] === missing
        A.data[chunk_key...] = fill(A.default, A.chunk_size...)
    end
    local_indices = ntuple(i -> mod(I[i] - 1, A.chunk_size[i]) + 1, N)
    A.data[chunk_key...][local_indices...] = value
end

function Base.iterate(A::ChunkedArray, state=1)
    if state > length(A)
        return nothing
    end
    I = ntuple(i -> div(state - 1, prod(A.dims[i+1:end])) % A.dims[i] + 1, length(A.dims))
    return (getindex(A, I...), state + 1)
end

function Base.show(io::IO, ::MIME"text/plain", a::ChunkedArray)
    print(io, join(a.dims, "x"))
    print(io, " ")
    print(io, typeof(a))
end