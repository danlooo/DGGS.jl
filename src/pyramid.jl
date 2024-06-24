
function DGGSPyramid(data::Dict{Int,DGGSLayer}, attrs=Dict{String,Any}())
    levels = data |> keys |> collect
    bands = data |> values |> first |> x -> x.bands
    dggs = data |> values |> first |> x -> x.dggs
    DGGSPyramid(data, attrs, levels, bands, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DGGSPyramid)
    printstyled(io, typeof(dggs); color=:white)
    println(io, "")
    println(io, "DGGS: $(dggs.dggs)")
    println(io, "Levels: $(dggs.levels)")

    println(io, "Non spatial axes:")
    for ax in dggs.data |> values |> first |> axes
        ax_name = DimensionalData.name(ax)
        startswith(String(ax_name), "q2di") && continue

        print(io, "  ")
        printstyled(io, ax_name; color=:red)
        print(io, " ")
        print(io, eltype(ax.val))
        println(io, "")
    end

    println(io, "Bands: ")
    for a in dggs.data |> first |> x -> values(x.second.data) |> collect
        print(io, "  ")
        Base.show(io, a)
        println(io, "")
    end
end

Base.getindex(dggs::DGGSPyramid, i::Integer) = dggs.data[i]

function open_dggs_pyramid(path::String)
    root_group = zopen(path)
    haskey(root_group.attrs, "_DGGS") || error("Zarr store is not in DGGS format")

    pyramid = Dict{Int,DGGSLayer}()
    for level in root_group.attrs["_DGGS"]["levels"]
        layer_ds = open_dataset(root_group.groups["$level"])
        pyramid[level] = DGGSLayer(layer_ds)
    end

    return DGGSPyramid(pyramid, root_group.attrs)
end


function write_dggs_pyramid(base_path::String, dggs::DGGSPyramid)
    mkdir(base_path)

    for level in dggs.levels
        level_path = "$base_path/$level"

        arrs = map((k, v) -> k => v.data, keys(dggs[level].data), values(dggs[level].data))

        attrs = dggs[level].attrs
        attrs["_DGGS"] = Dict{Symbol,Any}(
            key => getfield(dggs.dggs, key) for key in propertynames(dggs.dggs)
        )
        attrs["_DGGS"][:level] = level

        ds = Dataset(; properties=attrs, arrs...)
        savedataset(ds; path=level_path)
        JSON3.write("$level_path/.zgroup", Dict(:zarr_format => 2))
        Zarr.consolidate_metadata(level_path)

        # required for open_dggs_array using HTTP
        for band in keys(ds.cubes)
            attrs = JSON3.read("$level_path/$band/.zattrs", Dict{String,Any})
            attrs = merge(attrs, dggs[level][band].attrs)
            attrs["_DGGS"] = Dict{String,Any}(
                String(key) => getfield(dggs.dggs, key) for key in propertynames(dggs.dggs)
            )
            attrs["_DGGS"]["level"] = level
            JSON3.write("$level_path/$band/.zattrs", attrs)
            Zarr.consolidate_metadata("$level_path/$band")
        end
    end

    global_attrs = dggs.attrs
    global_attrs["_DGGS"] = Dict{String,Any}(String(key) => getfield(dggs.dggs, key) for key in propertynames(dggs.dggs))
    global_attrs["_DGGS"]["levels"] = dggs.levels

    JSON3.write("$base_path/.zattrs", global_attrs)
    JSON3.write("$base_path/.zgroup", Dict(:zarr_format => 2))
    Zarr.consolidate_metadata(base_path)
    return nothing
end

"position of last row or column in a quad matrix of that level"
width(level::Integer) = 2^(level - 1)

function aggregate_pentagon!(xout::AbstractArray, n::Integer, a::DGGSArray)
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
    res = map(i -> a.data.data[i...], children[n]) |> mean
    xout[1, 1] = res
end

function aggregate_hexagons!(xout::AbstractArray, xin::AbstractArray, n::Integer, a::DGGSArray)
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
        2 => first_row_fwd(3),
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
    padded_xin = hcat(xin.parent[col_paddings[n]...], xin) # i=1
    padded_xin = vcat(vcat([missing], xin.parent[row_paddings[n]...])', padded_xin) # j=1

    kernel = Float64[1 1 0; 1 2 1; 0 1 1] |> x -> x ./ sum(x)
    kernel_stride = 2
    offset_i = -1
    offset_j = -1
    padding = 1

    for j in 1:size(xout, 2)
        for i in 1:size(xout, 1)
            irange = (i-1)*kernel_stride+1+offset_i+padding:i*kernel_stride+1+offset_i+padding
            jrange = (j-1)*kernel_stride+1+offset_j+padding:j*kernel_stride+1+offset_j+padding

            data = view(padded_xin, irange, jrange)
            # kernel is already normalized, just sum instead of mean
            xout[i, j] = sum(data .* kernel)
        end
    end

    xout[1, 1] = missing  # pentagons are handled separateley (different kernel)
end

"""
Spatial hexagonal convolution in Q2DI index space matching levels of DGGRID ISEA4H grids
"""
function aggregate_dggs_layer(xout, xin, arr::DGGSArray)
    n = xin.indices[3]
    aggregate_hexagons!(xout, xin, n, arr)
    aggregate_pentagon!(xout, n, arr)
end

function to_dggs_pyramid(geo_ds::Dataset, level::Integer, args...; verbose=true, kwargs...)
    verbose && @info "Convert to DGGS layer"
    l = to_dggs_layer(geo_ds, level, args...; kwargs...)
    verbose && @info "Building pyramid"
    dggs = to_dggs_pyramid(l)
    return dggs
end

function to_dggs_pyramid(l::DGGSLayer)
    pyramid = Dict{Int,DGGSLayer}()
    pyramid[l.level] = l

    for coarser_level in l.level-1:-1:2
        finer_layer = pyramid[coarser_level+1]
        coarser_data = Dict{Symbol,DGGSArray}()
        for (k, arr) in finer_layer.data
            coarser_arr = mapCube(
                (xout, xin) -> aggregate_dggs_layer(xout, xin, arr),
                arr.data;
                indims=InDims(:q2di_i, :q2di_j),
                outdims=OutDims(
                    Dim{:q2di_i}(range(0; step=1, length=2^(coarser_level - 1))),
                    Dim{:q2di_j}(range(0; step=1, length=2^(coarser_level - 1))),
                    path=tempname() # disable inplace
                )
            )
            l_attrs_clean = filter(((k, v),) -> k != "_DGGS", l.attrs)
            attrs = merge(l_attrs_clean, arr.attrs)
            coarser_data[k] = DGGSArray(coarser_arr, attrs, k, coarser_level, finer_layer.dggs)
        end

        pyramid[coarser_level] = DGGSLayer(coarser_data, l.attrs)
    end

    return DGGSPyramid(pyramid, l.attrs)
end

to_dggs_pyramid(a::DGGSArray; kw...) = a |> DGGSLayer |> l -> to_dggs_pyramid(l; kw...)

function to_dggs_pyramid(raster::AbstractDimArray, level::Integer; kw...)
    arr = to_dggs_array(raster, level; kw...)
    dggs = to_dggs_pyramid(arr)
    return dggs
end