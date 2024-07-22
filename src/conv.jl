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

"""
clip range that may lie in padding to quad boundaries
"""
function clip(r::UnitRange, level)
    quad_size = 2^(level - 1)
    start = r.start <= 0 ? 1 : r.start
    stop = r.stop >= quad_size ? quad_size : r.stop
    return start:stop
end

function paddings()

end

function window(a::DGGSArray, i::Q2DI, radius::Int)
    radius >= 1 || error("range must stop with a positive number")
    radius == 1 && return a[i]

    quad_size = 2^(a.level - 1)
    window_size = 2 * radius - 1

    irange = i.i-(radius-1):i.i+(radius-1)
    jrange = i.j-(radius-1):i.j+(radius-1)

    # needs paddings from other quads?
    i_start_pad = i.i - radius < 0
    i_stop_pad = i.i + radius >= quad_size
    j_start_pad = i.j - radius < 0
    j_stop_pad = i.j + radius >= quad_size

    within_same_quad = all(.![i_start_pad, i_stop_pad, j_start_pad, j_stop_pad])
    pad_i_start = i_start_pad & !i_stop_pad & !j_start_pad & !j_stop_pad
    pad_j_start = !i_start_pad & !i_stop_pad & j_start_pad & !j_stop_pad

    if within_same_quad
        quad_range = i.i-(radius-1):i.i+(radius-1)
        res = a.data[q2di_n=i.n, q2di_i=quad_range, q2di_j=quad_range]
        return res
    end

    if pad_j_start
        padding_size = window_size - length(clip(jrange, a.level))

        # concat main and padding using MapCube
        # to preserve other dimensions e.g. time
        other_n = Dict(3 => 7, 4 => 8, 5 => 9, 6 => 10)
        jrange_padding = range(stop=quad_size, step=1, length=padding_size)
        padding = a.data[q2di_n=other_n[i.n], q2di_i=irange, q2di_j=jrange_padding]

        res = mapCube(
            # map over main data without padding
            a.data[
                q2di_n=i.n,
                q2di_i=clip(irange, a.level),
                q2di_j=clip(jrange, a.level)
            ],
            indims=InDims(:q2di_i, :q2di_j),
            outdims=OutDims(
                Dim{:q2di_i}(irange),
                Dim{:q2di_j}(jrange)
            )
        ) do xout, xin
            non_spatial_idx = xin.indices |> x -> x[3:length(x)]
            cur_padding = padding[:, :, non_spatial_idx...]
            xout .= hcat(cur_padding, xin)
        end
        return res
    end

    error("Not implemented")
end

"""
Return a center cell `i` and its neighbors within a distance of `disk_size`
"""
function disk(a::DGGSArray, i::Q2DI, disk_size::Integer; fill_value=NaN)
    mask = disk_mask(disk_size)
    w = window(a, i, disk_size)
    error("not implemented")
end

"""
Return cells having the same distance of `ring_size` of a given center cell `i`
"""
function ring(a::DGGSArray, i::Q2DI, ring_size::Integer; fill_value=NaN)
    mask = ring_mask(ring_size)
    w = window(a, i, ring_size)
    error("not implemented")
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