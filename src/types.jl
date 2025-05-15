struct Cell{T<:Integer}
    i::T
    j::T
    n::T
    resolution::Int8

    function Cell{T}(i, j, n, resolution) where {T<:Integer}
        0 <= n <= 4 || error("n=$n must be within 0 and 4")
        0 <= i <= 2 * 2^resolution - 1 || error("i=$i must be within 0 and $(2*2^resolution - 1 )")
        0 <= j <= 2^resolution - 1 || error("j=$j must be within 0 and $(2^resolution - 1)")

        new{T}(i, j, n, resolution)
    end

    function Cell(i, j, n, resolution)
        T = typeof(i)
        new{T}(T(i), T(j), T(n), resolution)
    end
end

struct DGGSArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me} <: AbstractDimArray{T,N,D,A}
    # DimArray fields
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me

    # DGGS fields
    resolution::Integer
    dggsrs::String

    function DGGSArray(
        data::A, dims::D, refdims::R, name::Na, metadata::Me,
        resolution::Integer, dggsrs::String
    ) where {D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me} where {T,N}
        dims_d = Dict([DD.name(x) => x for x in dims])

        :dggs_i in keys(dims_d) || error("Dimension :dggs_i must be present")
        :dggs_j in keys(dims_d) || error("Dimension :dggs_j must be present")
        :dggs_n in keys(dims_d) || error("Dimension :dggs_n must be present")

        map(x -> 0 <= x <= 2 * 2^resolution - 1, dims_d[:dggs_i]) |> all || error("Dimension dggs_i not in range")
        map(x -> 0 <= x <= 2^resolution - 1, dims_d[:dggs_j]) |> all || error("Dimension dggs_j not in range")
        map(x -> 0 <= x <= 4, dims_d[:dggs_n]) |> all || error("Dimension dggs_n not in range")

        new{T,N,D,R,A,Na,Me}(data, dims, refdims, name, metadata, resolution, dggsrs)
    end
end

"Set of DGGSArrays with aligned and shared dimensions at the same resolution."
struct DGGSDataset{K,T,N,L,D<:Tuple,R<:Tuple,LD,M,LM} <: AbstractDimStack{K,T,N,L}
    # DimStack fields
    data::L
    dims::D
    refdims::R
    layerdims::NamedTuple{K,LD}
    metadata::M
    layermetadata::NamedTuple{K,LM}

    # DGGS fields
    resolution::Integer

    function DGGSDataset(
        data, dims, refdims, layerdims::LD, metadata, layermetadata,
        resolution
    ) where LD<:NamedTuple{K} where K
        T = DD.data_eltype(data)
        N = length(dims)
        DGGSDataset{K,T,N}(
            data, dims, refdims, layerdims, metadata, layermetadata, resolution
        )
    end
    function DGGSDataset{K,T,N}(
        data::L, dims::D, refdims::R, layerdims::NamedTuple, metadata::M, layermetadata::NamedTuple, resolution
    ) where {K,T,N,L,D,R,M}
        new{K,T,N,L,D,R,typeof(values(layerdims)),M,typeof(values(layermetadata))}(data, dims, refdims, layerdims, metadata, layermetadata,
            resolution
        )
    end
end