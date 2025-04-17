function to_dggs_array(geo_array, resolution; agg_func::Function=mean, outtype=Float64, lon_name=:lon, lat_name=:lat, kwargs...)
    lon_dim = filter(x -> name(x) == lon_name, dims(geo_array))
    lat_dim = filter(x -> name(x) == lat_name, dims(geo_array))
    isempty(lon_dim) && error("Longitude dimension not found")
    isempty(lat_dim) && error("Latitude dimension not found")
    lon_dim = lon_dim[1]
    lat_dim = lat_dim[1]

    cells = [(lon, lat) for lon in lon_dim for lat in lat_dim] |>
            x -> to_cell(x, resolution) |>
                 x -> reshape(x, length(lon_dim), length(lat_dim))

    # get pixels to aggregate for each cell
    cell_coords = Dict{eltype(cells),Vector{CartesianIndex}}()
    for c in axes(cells, 2)
        for r in axes(cells, 1)
            cell = cells[r, c]
            if cell in keys(cell_coords)
                # multiple pixels per cell
                push!(cell_coords[cell], CartesianIndex(r, c))
            else
                cell_coords[cell] = [CartesianIndex(r, c)]
            end
        end
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
                res = agg_func([view(xin, c) for c in cell_coords] .|> x -> x[1])
                xout[cell.i+1, cell.j+1, cell.n+1] = res
            catch
                @warn "Unable to process cell" cell
            end
        end
    end

    return cell_array
end

function to_geo_array(dggs_array, lon_dim::DD.Dimension, lat_dim::DD.Dimension; kwargs...)
    resolution = dggs_array.dggs_j |> length |> log2 |> Int
    cells = [(lon, lat) for lon in lon_dim for lat in lat_dim] |>
            x -> to_cell(x, resolution) |>
                 x -> reshape(x, length(lon_dim), length(lat_dim))
    geo_array = mapCube(
        dggs_array,
        indims=InDims(
            dggs_array.dggs_i,
            dggs_array.dggs_j,
            dggs_array.dggs_n
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

DGGSArray(array::AbstractDimArray) = DGGSArray(array.data, dims(array), refdims(array), name(array), metadata(array))

Base.propertynames(::DGGSArray) = (:resolution, :dggsrs, fieldnames(DGGSArray)...)

function Base.getproperty(array::DGGSArray, v::Symbol)
    v == :resolution && return array.metadata["dggs_resolution"]
    v == :dggsrs && return array.metadata["dggs_dggsrs"]

    return getfield(array, v)
end

function Base.show(io::IO, ::MIME"text/plain", array::DGGSArray)
    println(io, "DGGSArray{$(typeof(array.data).name.name),$(eltype(array)),...} $(array.dggsrs) at resolution $(array.resolution)")
    println(io, "Additional dimensions:")
    for dim in array.dims
        name(dim) in [:dggs_i, :dggs_j, :dggs_n] && continue
        println(io, "   $(name(dim))")
    end
end