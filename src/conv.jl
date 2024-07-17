"""
Create a mask of cells arround a center cell within a radius `disk_size`
"""
function disk_mask(disk_size::Integer)
    mask_size = 2 * disk_size - 1
    res = [!(y - x > (disk_size - 1) || x - y > (disk_size - 1)) for x in 1:mask_size, y in 1:mask_size]
    return res
end

function ring_mask(ring_size::Integer)
    # a ring is the intersection of a disks of 2 consecutive disks
    mask_size = 2 * ring_size - 1
    a = disk_mask(ring_size)
    b = fill(false, mask_size, mask_size)
    b[2:mask_size-1, 2:mask_size-1] = disk_mask(ring_size - 1)
    res = Matrix{Bool}(a - b)
    return res
end

function is_in_same_quad(a::DGGSArray, i::Q2DI, disk_size::Integer)
    quad_size = 2^(a.level - 1)
    i_is_within = 1 <= i.i - (disk_size - 1) < i.i + (disk_size - 1) <= quad_size
    j_is_within = 1 <= i.j - (disk_size - 1) < i.j + (disk_size - 1) <= quad_size
    is_not_pentagon = !(i.i == i.j == 1)
    is_within = i_is_within && j_is_within && is_not_pentagon
    return is_within
end

function neighborhood(a::DGGSArray, mask::Matrix{Bool}, i::Q2DI, disk_size::Integer; fill_value=NaN)
    disk_size <= 0 && error("disk_size must be positive")
    disk_size == 1 && return a[i] |> collect
    is_in_same_quad(a, i, disk_size) || error("Multi quad spanning is not supported")
    range = i.i-(disk_size-1):i.i+(disk_size-1) # mask is a square matrix
    data = a.data[q2di_n=i.n, q2di_i=range, q2di_j=range] |> collect
    res = data .* map(x -> x ? x : fill_value, mask)
    return res
end

"""
Return a center cell `i` and its neighbors within a distance of `disk_size`
"""
function disk(a::DGGSArray, i::Q2DI, disk_size::Integer; fill_value=NaN)
    mask = disk_mask(disk_size)
    return neighborhood(a, mask, i, disk_size; fill_value=fill_value)
end

"""
Return cells having the same distance of `ring_size` of a given center cell `i`
"""
function ring(a::DGGSArray, i::Q2DI, ring_size::Integer; fill_value=NaN)
    mask = ring_mask(ring_size)
    return neighborhood(a, mask, i, ring_size; fill_value=fill_value)
end

"""
Spatial hexagonal convolution using an aggregation function
"""
function conv(f::Function, disk_size::Int, a::DGGSArray)
    error("unimplemented")
end

"""
Spatial hexagonal convolution using a kernel matrix
"""
function conv(kernel::Matrix, a::DGGSArray)
    error("unimplemented")
end

"""
Aggregating spatial hexagonal convolution used to build pyramids
"""
function conv(a::DGGSArray)
    error("unimplemented")
end