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
        new{T}(i, j, n, resolution)
    end
end

struct DGGSArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me} <: AbstractDimArray{T,N,D,A}
    # DimArray fields
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me

    function DGGSArray(
        data::A, dims::D, refdims::R, name::Na, metadata::Me
    ) where {D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me} where {T,N}
        dims_d = Dict([DD.name(x) => x for x in dims])

        :dggs_i in keys(dims_d) || error("Dimension :dggs_i must be present")
        :dggs_j in keys(dims_d) || error("Dimension :dggs_j must be present")
        :dggs_n in keys(dims_d) || error("Dimension :dggs_n must be present")

        resolution = dims_d[:dggs_i] |> length |> log2 |> Int
        dggsrs = "ISEA4D.P5"

        map(x -> 0 <= x <= 2 * 2^resolution - 1, dims_d[:dggs_i]) |> all || error("Dimension dggs_i not in range")
        map(x -> 0 <= x <= 2^resolution - 1, dims_d[:dggs_j]) |> all || error("Dimension dggs_j not in range")
        map(x -> 0 <= x <= 4, dims_d[:dggs_n]) |> all || error("Dimension dggs_n not in range")

        if metadata == DD.Dimensions.Lookups.NoMetadata
            new_metadata = Dict{String,Any}()
        elseif typeof(metadata) != Dict{String,Any}
            @info "Convert metadata names to String"
            new_metadata = Dict{String,Any}(metadata)
        end
        new_metadata["dggs_resolution"] = resolution
        new_metadata["dggs_dggsrs"] = dggsrs

        new{T,N,D,R,A,Na,Dict{String,Any}}(data, dims, refdims, name, new_metadata)
    end
end