using NearestNeighbors
using DataFrames
using Distances
using Tidier
using StaticArraysCore

abstract type AbstractIndex end

struct SeqNum <: AbstractIndex
    id::Integer
    level::Integer
end

struct HIndex <: AbstractIndex
    id::Int64
end

struct Geo{T<:Real} <: AbstractIndex
    lon::T
    lat::T
end

struct MemIndex <: AbstractIndex
    id::Int64
end

# see get_fragments_lookup
# keys: hindex_int range, value: number of unused cells before that fragment
# as a vector (one for each resolution)
const hindex_fragments_isea7h = Any[Dict{UnitRange{Int64},Int64}(0:11 => 0), Dict{UnitRange{Int64},Int64}(70:75 => 10, 0:5 => 0, 77:82 => 11, 49:54 => 7, 7:12 => 1, 56:61 => 8, 63:68 => 9, 21:26 => 3, 42:47 => 6, 28:33 => 4, 14:19 => 2, 35:40 => 5), Dict{UnitRange{Int64},Int64}(147:152 => 24, 7:41 => 1, 392:397 => 64, 105:139 => 17, 98:103 => 16, 0:5 => 0, 252:286 => 41, 497:531 => 81, 343:348 => 56, 294:299 => 48, 301:335 => 49, 49:54 => 8, 350:384 => 57, 399:433 => 65, 539:544 => 88, 203:237 => 33, 154:188 => 25, 546:580 => 89, 56:90 => 9, 196:201 => 32, 245:250 => 40, 441:446 => 72, 490:495 => 80, 448:482 => 73), Dict{UnitRange{Int64},Int64}(2058:2063 => 342, 3479:3723 => 578, 7:41 => 1, 1372:1377 => 228, 1036:1070 => 172, 1379:1413 => 229, 49:293 => 8, 1764:2008 => 293, 3430:3435 => 570, 2065:2099 => 343, 3780:3814 => 628, 0:5 => 0, 343:348 => 57, 686:691 => 114, 350:384 => 58, 1078:1322 => 179, 2401:2406 => 399, 3773:3778 => 627, 2793:3037 => 464, 2744:2749 => 456, 1715:1720 => 285, 2751:2785 => 457, 3822:4066 => 635, 3094:3128 => 514, 3437:3471 => 571, 2408:2442 => 400, 2450:2694 => 407, 1421:1665 => 236, 693:727 => 115, 1029:1034 => 171, 1722:1756 => 286, 392:636 => 65, 735:979 => 122, 2107:2351 => 350, 3087:3092 => 513, 3136:3380 => 521), Dict{UnitRange{Int64},Int64}(21952:23666 => 3657, 19257:19501 => 3208, 0:5 => 0, 9947:11661 => 1657, 2408:2442 => 401, 2450:2694 => 408, 14455:14699 => 2408, 7203:7208 => 1200, 21609:21614 => 3600, 12348:14062 => 2057, 16807:16812 => 2800, 49:293 => 8, 24353:26067 => 4057, 7210:7244 => 1201, 16814:16848 => 2801, 19208:19213 => 3200, 5145:6859 => 857, 17150:18864 => 2857, 12012:12046 => 2001, 7:41 => 1, 4851:5095 => 808, 4802:4807 => 800, 26411:26416 => 4400, 2401:2406 => 400, 14749:16463 => 2457, 9653:9897 => 1608, 24059:24303 => 4008, 9604:9609 => 1600, 21616:21650 => 3601, 7546:9260 => 1257, 14413:14447 => 2401, 4809:4843 => 801, 19551:21265 => 3257, 26418:26452 => 4401, 26754:28468 => 4457, 343:2057 => 57, 12005:12010 => 2000, 24010:24015 => 4000, 9611:9645 => 1601, 7252:7496 => 1208, 14406:14411 => 2400, 26460:26704 => 4408, 24017:24051 => 4001, 12054:12298 => 2008, 19215:19249 => 3201, 21658:21902 => 3608, 2744:4458 => 457, 16856:17100 => 2808), Dict{UnitRange{Int64},Int64}(120050:132054 => 20007, 36015:48019 => 6002, 33621:33655 => 5603, 134456:134461 => 22408, 151606:153320 => 25266, 151312:151556 => 25217, 0:5 => 0, 33614:33619 => 5602, 151263:151268 => 25209, 50421:50426 => 8403, 184877:184882 => 30811, 168077:168111 => 28011, 117992:119706 => 19664, 187278:199282 => 31211, 33957:35671 => 5659, 84378:86092 => 14062, 151270:151304 => 25210, 86436:98440 => 14405, 16807:16812 => 2801, 49:293 => 8, 168413:170127 => 28067, 84084:84328 => 14013, 153664:165668 => 25609, 67277:67521 => 11212, 16814:16848 => 2802, 184884:184918 => 30812, 100842:100847 => 16806, 84035:84040 => 14005, 17150:18864 => 2858, 50764:52478 => 8460, 100891:101135 => 16814, 168119:168363 => 28018, 7:41 => 1, 33663:33907 => 5610, 134505:134749 => 22416, 50428:50462 => 8404, 117698:117942 => 19615, 168070:168075 => 28010, 2401:14405 => 400, 84042:84076 => 14006, 100849:100883 => 16807, 136857:148861 => 22808, 134799:136513 => 22465, 52822:64826 => 8803, 343:2057 => 57, 184926:185170 => 30819, 67235:67269 => 11205, 19208:31212 => 3201, 117649:117654 => 19607, 134463:134497 => 22409, 69629:81633 => 11604, 67571:69285 => 11261, 170471:182475 => 28410, 101185:102899 => 16863, 16856:17100 => 2809, 67228:67233 => 11204, 50470:50714 => 8411, 103243:115247 => 17206, 117656:117690 => 19608, 185220:186934 => 30868)]

function normalize(x::Vector{<:Real}, digits=2)
    length(x) == digits && return x
    length(x) > digits && throw(ArgumentError("Need >= $(length(x)) digits for lossless normalization, but got $(digits)"))

    n_leading_zeros = Int8.(zeros(digits - length(x)))
    return vcat(n_leading_zeros, x)
end

function hindex_int(v::Vector)
    all(v .< 7) || throw(DomainError("All elements must be <7"))
    foldl((x, y) -> 7 * x + y, v) # conversion from base 10 to 7
end

function hindex_vector(i::Int, digits::Int=8)
    # conversion from base 7 to 10
    # each base 10 digit is one element of the final hierarchical position vector
    res = Vector{Int8}()
    while true
        d, r = divrem(i, 7)
        i = d
        push!(res, r)
        i == 0 && break
    end
    return res |> reverse |> x -> normalize(x, digits)
end

hindex_vector(i::HIndex) = hindex_vector(i.id)
hindex_int(i::HIndex) = i.id

function Base.show(io::IO, ::MIME"text/plain", i::HIndex)
    id_str = i.id |> DGGS.hindex_vector |> x -> join(x, ".") |> x -> lstrip(x, ['0', '.'])
    print(io, "HIndex $(i.id) $(id_str)")
end

function HIndex(v::Vector{Int64})
    v |> hindex_int |> HIndex
end

level(i::HIndex) = hindex_vector(i) |> join |> x -> lstrip(x, '0') |> length |> x -> x == 0 ? 1 : x

function MemIndex(i::HIndex)
    fragments = hindex_fragments_isea7h[level(i)]

    for (fragment, offset) in fragments
        if i.id in fragment
            return MemIndex(i.id - offset + 1)
        end
    end

    throw(ArgumentError("Invalid index"))
end