# adjecend functions e.g. to create constants that are not part of the DGGS module

using DGGS
using Distances

"""
Vector to store fragmentated index cohesiveley
e.g. index [1,2,3,7,8] to fragments 1:3 and 7:8
Missing dindices will be skiped to save memory in RAM and disk
"""
struct FragmentedVector{A<:Vector}
    data::A
    fragments
    offsets::Vector{Int64}
    length::Int64
end

function FragmentedVector(x::Array, y::Array)
    length(x) == length(y) || throw(ArgumentError("Lengths of x and y must be the same"))

    fragments = []
    global start = 1
    for pos in 2:length(x)
        if x[pos] > x[pos-1] + 1
            # new fragment detected
            push!(fragments, x[start]:x[pos-1])
            global start = pos
        end
    end
    push!(fragments, x[start]:x[length(x)])
    offsets = fragments[1:length(fragments)-1] |> x -> vcat(1, x) .|> length |> cumsum .|> x -> x - 1
    FragmentedVector(y, fragments, offsets, length(y))
end

FragmentedVector(y::Array) = FragmentedVector(y, [1:length(y)], [0], length(y)) # trivial case

Base.length(a::FragmentedVector) = a.length
Base.parent(a::FragmentedVector) = a.data
Base.eltype(a::FragmentedVector) = eltype(a.data)

function Base.getindex(a::FragmentedVector, i::Integer)
    for (f, fragment) in enumerate(a.fragments)
        if i in fragment
            return a.data[i-a.offsets[f]]
        end
    end
    throw(BoundsError(a, i))
end

function Base.setindex!(A::FragmentedVector, X, i::Integer)
    for fragment in a.fragments
        if i in fragment
            # position valid
            A.data[i] = X
            return
        end
    end
    throw(BoundsError(a, i))
end

function Base.show(io::IO, ::MIME"text/plain", a::FragmentedVector)
    println(io, "$(typeof(a))")
    print(io, "$(length(a.fragments)) fragments: ")
    println(io, a.fragments)
    print(io, "Data: ")
    print(io, a.data)
end



using StatsBase

function get_gaps(indices, level)
    index = indices[level]
    gaps = [0]
    gaps_count = 0
    for i in 1:maximum(index)
        if i in index
        else
            gaps_count += 1
        end
        push!(gaps, gaps_count)
    end
    return gaps
end

"""
Calculate cohesive fragments in Hindex of ISEA7H grid
Required to convert HIndex to MemIndex
"""
function get_fragments_lookup(levels=6)
    grids = DgGrid.(:hexagon, 7, :isea, 1:levels; distance=Haversine())
    indices = DGGS.seqnum2hindices(grids)

    res = []
    for (level, grid) in enumerate(grids)
        gaps = get_gaps(indices, level)

        x = sort(indices[level])
        y = repeat([1], length(grid))
        v = FragmentedVector(x, y)
        data = Dict(v.fragments[i] => gaps[v.fragments[i].start+1] for i = 1:length(v.fragments))
        push!(res, data)
    end
    res
end

fragments_lookup = get_fragments_lookup()
open("lookup.txt", "w") do f
    fragments_lookup |>
    repr |>
    x -> write(f, x)
end