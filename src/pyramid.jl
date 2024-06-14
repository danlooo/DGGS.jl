
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

function aggregate_dggs_layer(xout, xin, agg_func)
    fac = ceil(Int, size(xin, 1) / size(xout, 1))
    for j in axes(xout, 2)
        for i in axes(xout, 1)
            iview = ((i-1)*fac+1):min(size(xin, 1), (i * fac))
            jview = ((j-1)*fac+1):min(size(xin, 2), (j * fac))
            data = view(xin, iview, jview)
            xout[i, j] = agg_func(data)
        end
    end
end

function to_dggs_pyramid(geo_ds::Dataset, level::Integer, args...; verbose=true, kwargs...)
    verbose && @info "Convert to DGGS layer"
    l = to_dggs_layer(geo_ds, level, args...; kwargs...)
    verbose && @info "Building pyramid"
    dggs = to_dggs_pyramid(l)
    return dggs
end

function to_dggs_pyramid(l::DGGSLayer; agg_func::Function=filter_null(mean))
    pyramid = Dict{Int,DGGSLayer}()
    pyramid[l.level] = l

    for coarser_level in l.level-1:-1:2
        finer_layer = pyramid[coarser_level+1]
        coarser_data = Dict{Symbol,DGGSArray}()
        for (k, arr) in finer_layer.data
            coarser_arr = mapCube(
                (xout, xin) -> aggregate_dggs_layer(xout, xin, agg_func),
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