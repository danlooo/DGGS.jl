using DimensionalData.Dimensions: dimcolors

#
# grid topology
#

n_submats = 5
crs_geo = "EPSG:4326"
crs_isea = "+proj=isea +orient=isea +mode=plane +R=6371007.18091875"

# see https://www.youtube.com/watch?v=FWJOdMh8JQo&t=341s
θ_rotation = 60 * π / 180
rotation_matrix = [cos(θ_rotation) -sin(θ_rotation); sin(θ_rotation) cos(θ_rotation)]
θ_skew = 30 * π / 180
skew_matrix = [1 tan(θ_skew); 0 1]
trans_matrix = skew_matrix * rotation_matrix

trans = Proj.Transformation(crs_geo, crs_isea)
inv_trans = inv(trans)
epsilon = 1e-3 # move point from vertix into quad
# use icosahedron vertices to calculate the extent
x_rect_min, y_rect_min = trans_matrix * vcat(trans(58.2825256, -168.75)...)
x_rect_max, y_rect_max = trans_matrix * vcat(trans(0.0 + epsilon, -137.0325256 - epsilon)...)
rect_width = x_rect_max - x_rect_min
rect_height = y_rect_max - y_rect_min

#
# coordinate transformations
#

"""
Transform geographical coordinates (lat,lon) to cell ids (i,j,n)
Reverse operation of `to_geo`
"""
function to_cell(lon::Real, lat::Real, resolution)
    # sanity checks
    -180 <= lon <= 180 || error("lon must be within -180 and 180")
    -90 <= lat <= 90 || error("lat must be within -90 and 90")

    # edge cases
    lat == -90 && return Cell(2 * 2^resolution - 1, 0.5 * 2^resolution, 1, resolution)
    lat == 90 && return Cell(0, 0.5 * 2^resolution, 0, resolution)

    # project to ISEA
    trans = transformations[threadid()]
    x_isea, y_isea = trans(lat, lon)

    # rectify
    x_rect, y_rect = trans_matrix * [x_isea; y_isea]
    n_cell = (y_rect - y_rect_min) / rect_height * n_submats |> floor |> Int

    # affine transformation
    x_scaled = (x_rect - x_rect_min) / rect_width * (n_submats + 1)
    y_scaled = (y_rect - y_rect_min) / rect_height * n_submats

    x_offset = x_scaled - (n_cell + 1) + 1
    y_offset = y_scaled - n_cell

    # Discretize coordinates to cell ids
    i_cell = x_offset * 2^resolution |> floor |> Int
    j_cell = y_offset * 2^resolution |> floor |> Int

    # crop
    # TODO: review
    i_cell = clamp(i_cell, 0, 2 * 2^resolution - 1)
    j_cell = clamp(j_cell, 0, 2^resolution - 1)
    n_cell = clamp(n_cell, 0, 4)

    return Cell(i_cell, j_cell, n_cell, resolution)
end

"""
Multi-threaded version of to_cell
geo_points: Vector of (lon,lat) tuples
"""
function to_cell(geo_points::Vector{Tuple{A,B}}, resolution) where {A<:Real,B<:Real}
    # avoid overhead for small data
    length(geo_points) < 1e5 && return map(x -> to_cell(x..., resolution), geo_points)

    res = Vector{Cell}(undef, length(geo_points))
    Threads.@threads for i in eachindex(geo_points)
        lon, lat = geo_points[i][1], geo_points[i][2]
        cell = to_cell(lon, lat, resolution)
        res[i] = cell
    end

    return res
end

"""
Transform cell ids (i,j,n) to geographical coordinates (lat,lon)
Reverse of `to_cell_id`
"""
function to_geo(cell::Cell)
    # sanity checks are in Cell constructor

    # edge cases
    cell.i == 2 * 2^cell.resolution - 1 && return (0.0, -90.0)

    # Reverse Discretization
    x_scaled = (cell.i / 2^cell.resolution + cell.n)
    y_scaled = (cell.j / 2^cell.resolution + cell.n)

    # Reverse transformation
    x_rect = x_scaled / 6 * rect_width + x_rect_min
    y_rect = y_scaled / 5 * rect_height + y_rect_min

    x_isea, y_isea = inv(trans_matrix) * [x_rect; y_rect]

    # Reverse ISEA projection
    lat, lon = inv_transformations[threadid()](x_isea, y_isea)

    # Solve pole ambiguity
    # south pole already stable
    if lat >= 90
        return (0.0, 90.0)
    else
        return (lon, lat)
    end
end

function to_geo(cells::Vector{Cell{T}}) where {T<:Integer}
    # avoid overhead for small data
    length(cells) < 1e5 && return map(to_geo, cells)


    res = Vector{Tuple{Float64,Float64}}(undef, length(cells))
    Threads.@threads for i in eachindex(cells)
        geo_point = to_geo(cells[i])
        res[i] = geo_point
    end

    return res
end

to_geo(i, j, n, resolution) = Cell(i, j, n, resolution) |> to_geo

#
# Cell features
#

function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    print(io, typeof(cell), "(")
    printstyled(io, cell.i; color=dimcolors(1))
    print(io, ",")
    printstyled(io, cell.j; color=dimcolors(2))
    print(io, ",")
    printstyled(io, cell.n; color=dimcolors(3))
    print(io, ",")
    printstyled(io, cell.resolution; color=:white)
    print(io, ")")
end

function Base.getproperty(cell::Cell, name::Symbol)
    name == :id && return Int(cell)
    name == :geo && return to_geo(cell)
    name == :lon && return to_geo(cell)[1]
    name == :lat && return to_geo(cell)[2]

    return getfield(cell, name)
end

Base.propertynames(::Cell) = (:id, :geo, :lat, :lon, fieldnames(Cell)...)

Base.isless(a::Cell, b::Cell) = a.resolution == b.resolution ? Int(a) < Int(b) : error("Resolutions must be equal")

#
# conversions
#

"""
Int64 representation of the cell according to OGC

Bit Layout
==========
3 bits were added in each subsequent resolution to represent i and j coordinates.
This determines the cell ordering used for iteration and storage.
This linear index is designed to be compact, i.e. all cell numbers from 0 to length(cells)-1 were used.
In addition, the resolution is stored in the first 3 significant bits to be unambiguous across different spatial resolutions.
The traversing starts at the i coordinates longer side of the first of the 5 sub matrices, going line by line of the same sub matrix.
If all points within a sub matrix are visited, it seemingness continuous on the subsequent sub matrix. 

Most 3 significant bits: n ∈ [0,4]
Next resolution bits: j ∈ [0, 2^resolution-1]
Next resolution+1 bits: i ∈ [0, 2*2^resolution-1]
"""
function Base.Int(cell::Cell{T}) where {T<:Integer}
    res = cell.i + cell.j << (cell.resolution + 1) + cell.n << (2 * cell.resolution + 1)
    return res
end

function Cell(cell_int::Integer, resolution::Int)
    # needs resolution not stored in the integer index to ensure sequential id
    # resolution = Int((length(cell_bits) - 2) / 2) # does not work for i = 0 or j=0
    to_int(bits) = sum([bits[x] * 2^(x - 1) for x in length(bits):-1:1])

    bits = digits(cell_int, base=2, pad=2 * resolution + 4)
    n = bits[end-2:end] |> to_int
    j = bits[end-resolution-2:end-3] |> to_int
    i = bits[1:end-resolution-3] |> to_int
    cell = Cell(i, j, n, resolution)
    return cell
end