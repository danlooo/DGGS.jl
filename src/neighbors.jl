"""
Create a mask of cells arround a center cell within a radius `disk_size`
"""
function disk_mask(disk_size::Integer)
    mask_size = 2 * disk_size - 1
    res = [!(y - x > (disk_size - 1) || x - y > (disk_size - 1)) for x in 1:mask_size, y in 1:mask_size]
    return res
end

"""
Create a mask of cells having the same distance `ring_size` to a center cell
"""
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

function get_window_pad_nothing(a::DGGSArray, center::Q2DI, disk_size::Integer)
    # all cells of mask are within one quad
    irange = center.i-(disk_size-1):center.i+(disk_size-1)
    jrange = center.j-(disk_size-1):center.j+(disk_size-1)
    return a.data[q2di_n=center.n, q2di_i=irange, q2di_j=jrange]
end

function get_window_pad_j_start(a::DGGSArray, center::Q2DI, disk_size::Integer)
    irange = center.i-(disk_size-1):center.i+(disk_size-1)
    jrange = center.j-(disk_size-1):center.j+(disk_size-1)

    main = a.data[
        q2di_n=center.n,
        q2di_i=DGGS.clip(irange, a.level),
        q2di_j=DGGS.clip(jrange, a.level)
    ]
    non_spatial_axes = filter(x -> !startswith(String(DimensionalData.name(x)), "q2di"), a.data.axes)
    padding = YAXArray((
            main.q2di_i,
            Dim{:q2di_j}(jrange.start:0),
            non_spatial_axes...
        ),
        a.data[q2di_n=11, q2di_i=21:29, q2di_j=29:32],
        Dict()
    )
    padded = cat(padding, main, dims=Dim{:q2di_j}(jrange))
    return padded
end

function get_window_pad_i_start(a::DGGSArray, center::Q2DI, disk_size::Integer)
    irange = center.i-(disk_size-1):center.i+(disk_size-1)
    jrange = center.j-(disk_size-1):center.j+(disk_size-1)
    main = a.data[
        q2di_n=center.n,
        q2di_i=clip(irange, a.level),
        q2di_j=clip(jrange, a.level)
    ]
    non_spatial_axes = filter(x -> !startswith(String(DimensionalData.name(x)), "q2di"), a.data.axes)
    padding = YAXArray((
            Dim{:q2di_j}(9:17),
            Dim{:q2di_i}(-3:0),
            non_spatial_axes...
        ),
        a.data[q2di_n=6, q2di_i=reverse(16:24), q2di_j=29:32].data,
        Dict()
    )

    # permute dims to match main, needed for cat
    if length(non_spatial_axes) == 0
        padding = permutedims(padding, (2, 1))
    elseif length(non_spatial_axes) == 1
        padding = permutedims(padding, (2, 1, 3))
    else
        padding = permutedims(padding, (2, 1, 3:3+length(non_spatial_axes)...))
    end
    padded = cat(padding, main, dims=Dim{:q2di_i}(irange))
    return padded
end

function Base.getindex(a::DGGSArray, center::Q2DI, disk_size::Integer, type=:disk)
    type in [:disk, :ring, :window] || error("type not supported")

    mask = Dict(
        :disk => disk_mask(disk_size),
        :ring => ring_mask(disk_size),
        :window => (2 * disk_size - 1) |> x -> fill(true, x, x)
    )[type]

    irange = center.i-(disk_size-1):center.i+(disk_size-1)
    jrange = center.j-(disk_size-1):center.j+(disk_size-1)
    quad_size = width(a.level)
    i_is_in_same_quad = 1 <= irange.start <= irange.stop <= quad_size
    j_is_in_same_quad = 1 <= jrange.start <= jrange.stop <= quad_size

    if i_is_in_same_quad & j_is_in_same_quad
        window = get_window_pad_nothing(a, center, disk_size)
    elseif (irange.start < 1) & j_is_in_same_quad
        window = get_window_pad_i_start(a, center, disk_size)
        mask_size = size(mask)[1]
        mask = vcat(mask[1:disk_size-1, mask_size:-1:1], mask[disk_size:mask_size, :])
    elseif i_is_in_same_quad & (jrange.start < 1)
        window = get_window_pad_j_start(a, center, disk_size)
    else
        error("edge case not implemented")
    end

    masked = mapCube(
        (xout, xin) -> xout .= xin[mask],
        window,
        indims=InDims(:q2di_i, :q2di_j),
        outdims=OutDims(Dim{:q2di_k}(1:sum(mask)))
    )
    return masked
end

function Base.getindex(a::DGGSArray, lon::Real, lat::Real, disk_size::Integer, type=:disk)
    center = transform_points([(lon, lat)], a.level)[1]
    res = a[center, disk_size, type]
    return res
end