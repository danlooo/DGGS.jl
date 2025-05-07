function compute_cell_array(lon_dim, lat_dim, resolution)
    [(lon, lat) for lon in lon_dim, lat in lat_dim] |>
    x -> to_cell(x, resolution)
end

function to_dggs_array(geo_array, resolution; agg_func::Function=mean, outtype=Float64, lon_name=:X, lat_name=:Y, kwargs...)
    lon_dim = filter(x -> name(x) == lon_name, dims(geo_array))
    lat_dim = filter(x -> name(x) == lat_name, dims(geo_array))
    isempty(lon_dim) && error("Longitude dimension not found")
    isempty(lat_dim) && error("Latitude dimension not found")
    lon_dim = only(lon_dim)
    lat_dim = only(lat_dim)

    cells = compute_cell_array(lon_dim, lat_dim, resolution)

    # get pixels to aggregate for each cell
    cell_coords = Dict{eltype(cells),Vector{CartesianIndex{2}}}()
    for cI in CartesianIndices(cells)
        cell = cells[cI]
        current_cells = get!(() -> CartesianIndex{2}[], cell_coords, cell)
        push!(current_cells, cI)
    end

    # re-grid
    cell_array = mapCube(
        # mapCube can't find axes of other AbstractDimArrays e.g. Raster
        YAXArray(dims(geo_array), geo_array.data),
        indims=InDims(lon_dim, lat_dim),
        outdims=OutDims(
            Dim{:dggs_i}(0:(2*2^resolution-1)),
            Dim{:dggs_j}(0:(2^resolution-1)),
            Dim{:dggs_n}(0:4),
            outtype=outtype,
            kwargs...
        )) do xout, xin
        for (cell, cell_coords) in cell_coords
            try
                # view returns 0 dim array of pixels within the cell
                res = agg_func(view(xin, cell_coords))
                xout[cell.i+1, cell.j+1, cell.n+1] = res
            catch
                @warn "Unable to process cell" cell
            end
        end
    end

    return DGGSArray(cell_array, resolution, "ISEA4D.Penta")
end

function to_geo_array(dggs_array::DGGSArray, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    cells = compute_cell_array(lon_dim, lat_dim, dggs_array.resolution)
    geo_array = mapCube(
        dggs_array,
        indims=InDims(
            :dggs_i,
            :dggs_j,
            :dggs_n
        ),
        outdims=OutDims(lon_dim, lat_dim),
        kwargs...
    ) do xout, xin
        xout .= map(x -> xin[x.i+1, x.j+1, x.n+1], cells)
    end

    return geo_array
end

function to_geo_array(dggs_array, lon_range::AbstractRange, lat_range::AbstractRange; kwargs...)
    lon_dim = X(lon_range)
    lat_dim = Y(lat_range)
    to_geo_array(dggs_array, lon_dim, lat_dim; kwargs...)
end

#
# DGGSArray features
#

function DGGSArray(array::AbstractDimArray, resolution::Integer, dggsrs::String)
    return DGGSArray(
        array.data, dims(array), refdims(array), name(array), metadata(array),
        resolution, dggsrs
    )
end

function Base.show(io::IO, mime::MIME"text/plain", array::DGGSArray)
    println(io, "DGGSArray{$(eltype(array))} $(string(name(array)))")
    println(io, "DGGS: $(array.dggsrs) at resolution $(array.resolution)")

    if length(array.dims) > 3
        println(io, "Additional dimensions:")
        for dim in array.dims
            name(dim) in [:dggs_i, :dggs_j, :dggs_n] && continue
            print(io, "   ")
            DD.Dimensions.print_dimname(io, dim)
            print(io, " $(minimum(dim):step(dim):maximum(dim))")
        end
        println(io, "")
    else
        println(io, "No additional dimensions")
    end

    if length(array.metadata) > 0
        println(io, "Meta data:")
        for (key, value) in array.metadata
            println(io, "   $key: $value")
        end
    else
        println(io, "No meta data")
    end
end

"rebuild immutable objects with new field values. Part of any AbstractDimArray."
function DD.rebuild(
    dggs_array::DGGSArray, data::AbstractArray, dims::Tuple, refdims::Tuple, name, metadata
)
    DGGSArray(data, dims, refdims, name, metadata, dggs_array.resolution, dggs_array.dggsrs)
end


#
# Subset array
#

function Base.getindex(dggs_array::DGGSArray, cell::Cell)
    cell.resolution == dggs_array.resolution || error("Resolutions of cell and array must be the same")
    # test on DiskArrays
    return Base.getindex(DimArray(dggs_array), dggs_i=cell.i + 1, dggs_j=cell.j + 1, dggs_n=cell.n + 1)
end
