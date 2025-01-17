
function DGGSPyramid(data::AbstractDict{Int,DGGSLayer}, attrs=Dict{String,Any}())
    levels = data |> keys |> collect
    dggs = data |> values |> first |> x -> x.dggs
    DGGSPyramid(data, attrs, levels, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DGGSPyramid)
    printstyled(io, typeof(dggs); color=:white)
    println(io, "")
    println(io, "DGGS: $(dggs.dggs)")
    println(io, "Levels: $([x for x in dggs.levels])")

    show_axes(io, dggs.data |> values |> first |> axes)
    show_arrays(io, dggs.data |> first |> x -> values(x.second.data) |> collect)
end

Base.getindex(dggs::DGGSPyramid, level::Integer) = dggs.data[level]

function Base.getindex(p::DGGSPyramid; kwargs...)
    level = get(kwargs, :level, nothing)
    isnothing(level) && error("level not provided")
    length(kwargs) == 1 && return p[level]

    id = get(kwargs, :id, nothing)
    isnothing(id) && error("Array id not provided.")

    center = get(kwargs, :center, nothing)
    lon = get(kwargs, :lon, nothing)
    lat = get(kwargs, :lat, nothing)
    radii = get(kwargs, :radii, nothing)

    args = filter(!isnothing, (center, lon, lat, radii))
    kwargs = filter(x -> !(x.first in [:id, :level, :center, :lat, :lon, :radii]), kwargs)

    if isempty(args) & isempty(kwargs)
        return p.data[level][id]
    else
        return getindex(p.data[level][id], args...; kwargs...)
    end
end

function open_dggs_pyramid(path::String)
    root_group = zopen(path)
    haskey(root_group.attrs, "dggs_id") || error("Zarr store is not in DGGS format")

    pyramid = Dict{Int,DGGSLayer}()
    for level in sort(root_group.attrs["dggs_levels"])
        layer_ds = open_dataset(root_group.groups["$level"])
        pyramid[level] = DGGSLayer(layer_ds)
    end
    pyramid = sort!(OrderedDict(pyramid))

    return DGGSPyramid(pyramid, root_group.attrs)
end


function write_dggs_pyramid(base_path::String, dggs::DGGSPyramid)
    mkdir(base_path)

    for level in dggs.levels
        level_path = "$base_path/$level"

        arrs = map((k, v) -> k => v.data, keys(dggs[level].data), values(dggs[level].data))

        attrs = dggs[level].attrs
        attrs["dggs_level"] = level

        ds = Dataset(; properties=attrs, arrs...)
        savedataset(ds; path=level_path)
        JSON3.write("$level_path/.zgroup", Dict(:zarr_format => 2))
        Zarr.consolidate_metadata(level_path)

        # required for open_dggs_array using HTTP
        for band in keys(ds.cubes)
            attrs = JSON3.read("$level_path/$band/.zattrs", Dict{String,Any})
            attrs = merge(attrs, dggs[level][band].attrs)
            attrs["dggs_level"] = level
            JSON3.write("$level_path/$band/.zattrs", attrs)
            Zarr.consolidate_metadata("$level_path/$band")
        end
    end

    global_attrs = dggs.attrs
    global_attrs["dggs_levels"] = dggs.levels

    JSON3.write("$base_path/.zattrs", global_attrs)
    JSON3.write("$base_path/.zgroup", Dict(:zarr_format => 2))
    Zarr.consolidate_metadata(base_path)
    return nothing
end

"position of last row or column in a quad matrix of that level"
width(level::Integer) = 2^(level - 1)

function aggregate_pentagon!(xout::AbstractArray, xin::AbstractArray, n::Integer, a::DGGSArray; agg_type::Symbol=:round, lck=ReentrantLock())
    m = width(a.level)

    # position of children in Q2DI space i.e. (i,j,n)
    # see ../docs/src/assets/pentagon-children-q2di.png
    first_square(n) = [(1, 1, n), (1, 2, n), (2, 1, n), (2, 2, n)]
    children = Dict(
        1 => [(1, 1, 1), (1, m, 2), (1, m, 3), (1, m, 4), (1, m, 5), (1, m, 6)],
        2 => vcat(first_square(2), [(m, m, 6), (1, m, 11)]),
        3 => vcat(first_square(3), [(m, m, 2), (1, m, 7)]),
        4 => vcat(first_square(4), [(m, m, 3), (1, m, 8)]),
        5 => vcat(first_square(5), [(m, m, 4), (1, m, 9)]),
        6 => vcat(first_square(6), [(m, m, 5), (1, m, 10)]),
        7 => vcat(first_square(7), [(m, 1, 2), (m, m, 11)]),
        8 => vcat(first_square(8), [(m, m, 7), (m, 1, 3)]),
        9 => vcat(first_square(9), [(m, m, 8), (m, 1, 4)]),
        10 => vcat(first_square(10), [(m, m, 9), (m, 1, 5)]),
        11 => vcat(first_square(11), [(m, m, 10), (m, 1, 6)]),
        12 => [(1, 1, 12), (m, 1, 11), (m, 1, 10), (m, 1, 9), (m, 1, 8), (m, 1, 7)]
    )

    vals = map(children[n]) do (i, j, n)
        # select current slice at given spatial point
        # spatial dims are always at the beginning of idx
        idx = collect(xin.indices)
        idx[1] = i
        idx[2] = j
        idx[3] = n
        @lock lck a.data.data[idx...]
    end
    res = filter_null(mean)(vals) |> Dict(:round => round, :convert => identity, :identity => identity)[agg_type]
    xout[1, 1] = res
end

function aggregate_hexagons!(xout::AbstractArray, xin::AbstractArray, n::Integer, a::DGGSArray; agg_type::Symbol=:round, lck=ReentrantLock())
    (1 <= n <= 12) || error("Quad number n must be between 1 and 12")
    n in [1, 12] && return # first and last quad only contain one pentagon

    # padding to fill cell of other quads
    # Then, the convolution can run on each quad independently

    m = width(a.level)
    first_row_fwd(n) = (1, :, n)
    first_row_rev(n) = (1, m:-1:1, n)
    last_row_fwd(n) = (m, :, n)
    last_row_rev(n) = (m, m:-1:1, n)

    first_col_fwd(n) = (:, 1, n)
    first_col_rev(n) = (m:-1:1, 1, n)
    last_col_fwd(n) = (:, m, n)
    last_col_rev(n) = (m:-1:1, m, n)

    col_paddings = Dict(
        2 => last_row_rev(6),
        3 => last_col_fwd(7),
        4 => last_col_fwd(8),
        5 => last_col_fwd(9),
        6 => last_col_fwd(10),
        7 => last_row_rev(11),
        8 => last_row_rev(7),
        9 => last_row_rev(8),
        10 => last_row_rev(9),
        11 => last_row_rev(10),
    )
    row_paddings = Dict(
        2 => last_col_rev(6),
        3 => last_col_rev(2),
        4 => last_col_rev(3),
        5 => last_col_rev(4),
        6 => last_col_rev(5),
        7 => last_row_fwd(2),
        8 => last_row_fwd(3),
        9 => last_row_fwd(4),
        10 => last_row_fwd(5),
        11 => last_row_fwd(6),
    )

    col_idx = collect(xin.indices)
    col_idx[1] = col_paddings[n][1]
    col_idx[2] = col_paddings[n][2]
    col_idx[3] = col_paddings[n][3]

    row_idx = collect(xin.indices)
    row_idx[1] = row_paddings[n][1]
    row_idx[2] = row_paddings[n][2]
    row_idx[3] = row_paddings[n][3]

    padded_xin = @lock lck hcat(a.data[col_idx...], xin)
    padded_xin = @lock lck vcat(vcat([missing], a.data[row_idx...])', padded_xin)

    kernel = Float64[1 1 0; 1 2 1; 0 1 1] |> x -> x ./ sum(x)
    kernel_stride = 2
    offset_i = -1
    offset_j = -1
    padding = 1

    for j in axes(xout, 2)
        for i in axes(xout, 1)
            irange = (i-1)*kernel_stride+1+offset_i+padding:i*kernel_stride+1+offset_i+padding
            jrange = (j-1)*kernel_stride+1+offset_j+padding:j*kernel_stride+1+offset_j+padding

            data = view(padded_xin, irange, jrange)
            if all(ismissing.(data) .| isnan.(data))
                xout[i, j] = missing
                continue
            end

            # kernel is already normalized, just sum instead of mean
            xout[i, j] = filter_null(sum)(data .* kernel) |> Dict(:round => round, :convert => identity, :identity => identity)[agg_type]
        end
    end

    xout[1, 1] = missing  # pentagons are handled separateley (different kernel)
end

function to_dggs_pyramid(geo_ds::Dataset, level::Integer, args...; verbose=true, base_path=tempname(), agg_type=:round, kwargs...)
    verbose && @info "Convert to DGGS layer"
    l = to_dggs_layer(geo_ds, level, args...; agg_type=agg_type, kwargs...)
    verbose && @info "Building pyramid"
    dggs = to_dggs_pyramid(l; base_path=base_path, agg_type=agg_type)
    return dggs
end

function to_dggs_pyramid(l::DGGSLayer; base_path=tempname(), agg_type::Symbol=:round)
    pyramid = Dict{Int,DGGSLayer}()
    pyramid[l.level] = l

    for coarser_level in l.level-1:-1:2
        finer_layer = pyramid[coarser_level+1]
        coarser_data = Dict{Symbol,DGGSArray}()
        for (k, arr) in finer_layer.data
            if agg_type == :round && all(Base.uniontypes(eltype(arr.data)) .<: Union{AbstractFloat,Missing})
                # no rounding needed
                agg_type = :identity
            end

            # Some arrays e.g. NetDF4 over HDF5 are not thread-safe
            lck = ReentrantLock()
            coarser_arr = mapCube(
                arr.data;
                indims=InDims(:q2di_i, :q2di_j),
                outdims=OutDims(
                    Dim{:q2di_i}(range(0; step=1, length=2^(coarser_level - 1))),
                    Dim{:q2di_j}(range(0; step=1, length=2^(coarser_level - 1))),
                    path=joinpath(base_path, "$(coarser_level)/$(k)")
                )
            ) do xout, xin
                n = xin.indices[3]
                aggregate_hexagons!(xout, xin, n, arr; agg_type=agg_type, lck=lck)
                aggregate_pentagon!(xout, xin, n, arr; agg_type=agg_type, lck=lck)
            end
            attrs = deepcopy(arr.attrs)
            attrs["dggs_level"] = coarser_level
            coarser_arr = YAXArray(coarser_arr.axes, coarser_arr.data, attrs)
            coarser_data[k] = DGGSArray(coarser_arr, attrs, k, coarser_level, finer_layer.dggs)
        end

        pyramid[coarser_level] = DGGSLayer(coarser_data, l.attrs)
    end
    p = sort!(OrderedDict(pyramid))
    @warn "Pyramid may not be persistent on disk. Use write_dggs_pyramid for this."
    return DGGSPyramid(p, l.attrs)
end

to_dggs_pyramid(a::DGGSArray; kw...) = a |> DGGSLayer |> l -> to_dggs_pyramid(l; kw...)

function to_dggs_pyramid(raster::AbstractDimArray, level::Integer; agg_type=:round, kwargs...)
    arr = to_dggs_array(raster, level; agg_type=agg_type, kwargs...)
    dggs = to_dggs_pyramid(arr; agg_type=agg_type)
    return dggs
end