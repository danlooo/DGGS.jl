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

function DGGSArrayPyramid(cell_cube::DGGSArray; agg_func=filter_null(mean))
    pyramid = Dict{Int,DGGSArray}()
    pyramid[cell_cube.level] = cell_cube

    for coarser_level in cell_cube.level-1:-1:2
        coarser_cell_array = mapCube(
            (xout, xin) -> aggregate_cell_cube(xout, xin; agg_func=agg_func),
            pyramid[coarser_level+1].data,
            indims=InDims(:q2di_i, :q2di_j),
            outdims=OutDims(
                Dim{:q2di_i}(range(0; step=1, length=2^(coarser_level - 1))),
                Dim{:q2di_j}(range(0; step=1, length=2^(coarser_level - 1)))
            )
        )
        coarser_cell_array = YAXArray(coarser_cell_array.axes, coarser_cell_array.data, cell_cube.data.properties)
        coarser_cell_cube = DGGSArray(coarser_cell_array, coarser_level)
        pyramid[coarser_level] = coarser_cell_cube
    end

    return DGGSArrayPyramid(pyramid)
end

function DGGSArrayPyramid(data::AbstractArray{<:Number}, lon_range::AbstractVector, lat_range::AbstractVector, level::Integer; kwargs...)
    raster = DimArray(data, (X(lon_range), Y(lat_range)))
    cell_cube = to_dggs_array(raster, level; kwargs...)
    DGGSArrayPyramid(cell_cube)
end

function DGGSArrayPyramid(data::AbstractArray{<:Number}, lon_range::DimensionalData.XDim, lat_range::DimensionalData.YDim, level::Integer; kwargs...)
    raster = DimArray(data, (lon_range, lat_range))
    cell_cube = to_dggs_array(raster, level; kwargs...)
    DGGSArrayPyramid(cell_cube)
end

function GridSystem_url(url)
    url = replace(url, r"/+$" => "")
    tmpfile = download("$url/.zattrs")
    str = read(tmpfile, String)
    rm(tmpfile)
    attr_dict = JSON3.read(str)
    levels = attr_dict.grid.resolutions.spatial.levels
    length(levels) > 0 || error("No resolution level detected")

    pyramid = Dict{Int,DGGSArray}()
    for level in levels
        cell_array = open_dataset("$url/$level") |> Cube
        pyramid[level] = DGGSArray(cell_array, level)
    end

    return DGGSArrayPyramid(pyramid)
end

function GridSystem_local(path::String)
    isdir(path) || error("path '$path' must be a directory")

    attr_dict = JSON3.read("$path/.zattrs")
    levels = attr_dict.grid.resolutions.spatial.levels
    length(levels) > 0 || error("No resolution level detected")

    pyramid = Dict{Int,DGGSArray}()
    for level in levels
        cell_array = Cube("$path/$level")
        pyramid[level] = DGGSArray(cell_array, level)
    end

    return DGGSArrayPyramid(pyramid)
end

function DGGSArrayPyramid(path::String)
    startswith(path, r"http[s]?://") && return GridSystem_url(path)
    isdir(path) && return GridSystem_local(path)

    error("Path must be either a valid directory path or an URL")
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DGGSArrayPyramid)
    println(io, "DGGS Pyramid")
    println(io, "Levels: $(join(dggs.data |> keys |> collect |> sort, ","))")
    Base.show(io, "text/plain", dggs.data |> values |> first |> x -> x.data.axes)
end

function to_dggs_dataset_pyramid(geo_ds::Dataset, max_level::Int, cell_ids::DimArray{Q2DI{Int64},2})
    array_pyramids = OrderedDict{Symbol,DGGSArrayPyramid}()
    # build dggs for each variable
    # to ensure pyramid building in DGGS space
    for (var, geo_cube) in geo_ds.cubes
        axs = []
        for ax in geo_cube.axes
            if ax isa Dim{:lon}
                ax = X(ax.val .- 180)
                push!(axs, ax)
            elseif ax isa Dim{:lat}
                ax = Y(ax.val |> collect)
                push!(axs, ax)
            else
                push!(axs, ax)
            end
        end
        geo_cube = YAXArray(Tuple(axs), geo_cube.data, geo_cube.properties)
        cell_cube = to_dggs_array(geo_cube, max_level; cell_ids=cell_ids)
        array_pyramids[var] = cell_cube |> DGGSArrayPyramid
    end

    properties = geo_ds.properties
    properties["_DGGS"] = Dict(
        "index" => "Q2DI",
        "aperture" => 4,
        "rotation_lon" => 11.25,
        "polyhedron" => "icosahedron",
        "name" => "DGGRID ISEA4H Q2DI",
        "radius" => 6371007.180918475,
        "polygon" => "hexagon",
        "rotation_lat" => 58.2825,
        "projection" => "+isea",
        "rotation_azimuth" => 0
    )

    dggs = Dict{Integer,DGGSDataset}()
    for level in 2:max_level
        cubes = Dict{Symbol,YAXArray}()
        for (k, ar) in array_pyramids
            cubes[k] = array_pyramids[k][level].data
        end
        properties["_DGGS"]["level"] = level
        ds = Dataset(; properties=properties, cubes...) |> DGGSDataset
        dggs[level] = ds
    end
    dggs = DGGSDatasetPyramid(dggs)
    return dggs
end

function Base.write(base_path, dggs::DGGSDatasetPyramid)
    @infiltrate
end

function open_dggs_dataset(path::String)
    root_group = zopen(path)
    dggs_datasets = Dict()
    for (k, zarr_group) in sort(root_group.groups)
        ds = open_dataset(zarr_group)
        level = ds.properties["_DGGS"]["level"]
        dggs_datasets[level] = ds |> DGGSDataset
    end
    dggs = dggs_datasets |> DGGSDatasetPyramid
    return dggs
end

Base.getindex(dggs::DGGSDatasetPyramid, level::Int) = dggs.data[level]
Base.getindex(dggs::DGGSArrayPyramid, level::Int) = dggs.data[level]
Base.setindex!(dggs::DGGSArrayPyramid, cell_cube::DGGSArray, level::Int) = dggs.data[level] = cell_cube

DGGSArray(dggs::DGGSArrayPyramid) = dggs[dggs.data|>keys|>maximum]
Makie.plot(dggs::DGGSArrayPyramid, args...; kw...) = Makie.plot(DGGSArray(dggs), args...; kw...)

function get_axes(dggs::DGGSDatasetPyramid)
    level = dggs.data[dggs.data|>keys|>first]
    return level.data.axes |> values
end

function Base.show(io::IO, ::MIME"text/plain", dggs::DGGSDatasetPyramid)
    level = dggs.data[dggs.data|>keys|>first]
    axes = get_axes(dggs)

    println(io, typeof(dggs))

    print(io, "DGGS: ")
    printstyled(io, dggs.grid.name; color=:white)
    println(io)

    print(io, "Levels: ")
    printstyled(io, join(dggs.data |> keys |> collect |> sort, ", "); color=:white)
    println(io)

    println(io, "Axes:")
    for ax in axes
        startswith(String(name(ax)), "q2di_") && continue
        print(io, "  ")
        Base.show(io, "text/plain", ax)
        println(io)
    end

    println(io, "Variables:")
    for (k, arr) in level.data.cubes
        print(io, "  ")
        get(arr.properties, "standard_name", k) |> x -> printstyled(io, x; color=:blue)
        get(arr.properties, "units", "units undefined") |> x -> printstyled(io, " $x"; color=:white)
        print(io, " ")
        print(io, setdiff(name(arr.axes), (:q2di_n, :q2di_i, :q2di_j)) |> Tuple)
        print(io, " $(eltype(arr.data))")
        println(io)
    end

    println(io, "$(length(dggs.attrs)) Attributes")
end

function Base.getproperty(dggs::DGGSDatasetPyramid, v::Symbol)
    if v == :attrs
        return dggs.data |> values |> first |> x -> x.data.properties
    elseif v == :grid
        return dggs.attrs["_DGGS"] |> DGGSGridSystem
    elseif v == :axes
        return get_axes(dggs)
    else
        return getfield(dggs, v)
    end
end

Base.propertynames(::DGGSDatasetPyramid) = (:attrs, :data, :grid, :axes)
