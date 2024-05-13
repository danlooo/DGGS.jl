
function DGGSPyramid(data::Dict{Int,DGGSLayer}, attrs=Dict{String,Any}())
    levels = data |> keys |> collect
    bands = data |> values |> first |> x -> x.bands
    dggs = data |> values |> first |> x -> x.dggs
    DGGSPyramid(data, attrs, levels, bands, dggs)
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DGGSPyramid)
    println(io, "$(typeof(dggs))")
    println(io, "DGGS: $(dggs.dggs)")
    println(io, "Levels: $(dggs.levels)")
    println(io, "Bands: $(dggs.bands)")
end

Base.getindex(dggs::DGGSPyramid, i::Integer) = dggs.data[i]


function open_pyramid(path::String)
    root_group = zopen(path)
    haskey(root_group.attrs, "_DGGS") || error("Zarr store is not in DGGS format")

    pyramid = Dict{Int,DGGSLayer}()
    for level in root_group.attrs["_DGGS"]["levels"]
        layer_ds = open_dataset(root_group.groups["$level"])
        pyramid[level] = DGGSLayer(layer_ds)
    end

    return DGGSPyramid(pyramid, root_group.attrs)
end


function write_pyramid(base_path::String, dggs::DGGSPyramid)
    mkdir(base_path)

    for level in dggs.levels
        level_path = "$base_path/$level"
        arrs = map((k, v) -> k => v.data, keys(dggs[level].data), values(dggs[level].data))
        ds = Dataset(; properties=dggs[level].attrs, arrs...)
        savedataset(ds; path=level_path)
        JSON3.write("$level_path/.zgroup", Dict(:zarr_format => 2))
        Zarr.consolidate_metadata(level_path)

        # required for open_array using HTTP
        for band in keys(ds.cubes)
            Zarr.consolidate_metadata("$level_path/$band")
        end
    end

    JSON3.write("$base_path/.zattrs", dggs.attrs)
    JSON3.write("$base_path/.zgroup", Dict(:zarr_format => 2))
    Zarr.consolidate_metadata(base_path)
    return
end

function aggregate_cell_cube(xout, xin; agg_func=filter_null(mean))
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


function to_pyramid(raster::AbstractDimArray, level::Integer)
    data = Dict{Int,DGGSLayer}

    arr = to_array(raster, level)
    # TODO: to_array |> to_pyramid
    error("Not implemented")
end