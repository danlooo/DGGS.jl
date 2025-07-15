module DGGSMakie

using DGGS
using Makie
using DimensionalData
import DimensionalData as DD
using Makie.GeometryBasics
using Dates
using Extents
using Infiltrator

function get_resolution(dggs_pyramid::DGGSPyramid, lon_dim::X, lat_dim::Y)
    dggs_pyramid.dggsrs == "ISEA4D.Penta" || error("DGGSRS $(dggs_pyramid.dggsr) not implemented")

    global_lon_length = length(-180:Float64(lon_dim.val.step):180)
    global_lat_length = length(-90:Float64(lat_dim.val.step):90)
    global_pixels = global_lon_length * global_lat_length
    optimal_resolution = sqrt(global_pixels / 10) |> log2 |> ceil |> Int
    available_resolution = minimum([optimal_resolution, DGGS.get_resolutions(dggs_pyramid) |> maximum])
    return available_resolution
end

function get_texture(dggs_pyramid::DGGSPyramid, lon_dim::X, lat_dim::Y, layer::Symbol)
    resolution = get_resolution(dggs_pyramid, lon_dim, lat_dim)
    dggs_array = getproperty(dggs_pyramid[resolution], layer)
    return get_texture(dggs_array, lon_dim, lat_dim)
end

function get_texture(dggs_pyramid::DGGSPyramid, lon_dim::X, lat_dim::Y,
    red_layer::Symbol, green_layer::Symbol, blue_layer::Symbol,
    scale_factor::Real=1, offset::Real=0
)
    resolution = get_resolution(dggs_pyramid, lon_dim, lat_dim)
    dggs_ds = dggs_pyramid[resolution]
    return get_texture(
        dggs_ds, lon_dim, lat_dim,
        red_layer, green_layer, blue_layer,
        scale_factor, offset
    )
end

function get_texture(dggs_ds::DGGSDataset, lon_dim::X, lat_dim::Y, layer::Symbol)
    dggs_array = getproperty(dggs_ds, layer)
    return get_texture(dggs_array, lon_dim, lat_dim)
end

function get_texture(
    ds::DGGSDataset, lon_dim::X, lat_dim::Y,
    red_layer::Symbol, green_layer::Symbol, blue_layer::Symbol,
    scale_factor::Real=1, offset::Real=0
)
    ds_rgb = DGGSDataset(getproperty(ds, red_layer), getproperty(ds, green_layer), getproperty(ds, blue_layer))
    geo_ds = to_geo_dataset(ds_rgb, lon_dim, lat_dim)

    texture = map(CartesianIndices((lon_dim, lat_dim))) do i
        r = getproperty(geo_ds, red_layer)[i] * scale_factor + offset
        g = getproperty(geo_ds, green_layer)[i] * scale_factor + offset
        b = getproperty(geo_ds, blue_layer)[i] * scale_factor + offset

        if ismissing(r) || ismissing(g) || ismissing(b)
            Makie.RGBf(1, 1, 1)
        else
            Makie.RGBf(r, g, b)
        end
    end
    texture
end

function get_texture(dggs_array::DGGSArray, lon_dim::X, lat_dim::Y)
    return to_geo_array(dggs_array, lon_dim, lat_dim)
end

function Makie.plot(
    dggs::Union{DGGSArray,DGGSDataset,DGGSPyramid},
    args...
    ;
    extent=dggs.bbox,
    resolution_scale::Real=1
)
    fig = Figure()
    ax = Axis(fig[1, 1], limits=(extent.X..., extent.Y...), xlabel="longitude [°]", ylabel="latitude [°]")
    ax.aspect = DataAspect()

    data = @lift begin
        # get limits observable once. It might change during the update
        lims = $(ax.finallimits)
        lon_min, lat_min = lims.origin
        lon_max, lat_max = lims.origin .+ lims.widths

        lon_min = clamp(lon_min, -180, 180)
        lon_max = clamp(lon_max, -180, 180)
        lat_min = clamp(lat_min, -90, 90)
        lat_max = clamp(lat_max, -90, 90)

        lon_length, lat_length = $(fig.scene.viewport).widths
        lon_range = range(lon_min, lon_max, length=lon_length * resolution_scale)
        lat_range = range(lat_min, lat_max, length=lat_length * resolution_scale)
        lon_dim = X(lon_range)
        lat_dim = Y(lat_range)

        get_texture(dggs, lon_dim::X, lat_dim::Y, args...)
    end

    last_update_limits = Observable(ax.finallimits[])
    last_update_viewport_widths = Observable(fig.scene.viewport[].widths)

    # can't plot on every limit change
    # instead, plot one frame after the other if limits changed
    @async while true
        # skip update if limits did not change
        if ax.finallimits[] == last_update_limits[] && fig.scene.viewport[].widths == last_update_viewport_widths[]
            # FPS limit to waste less computation
            sleep(Millisecond(30))
            continue
        end

        # delay plotting if small zoom/pan detected
        change_frac = maximum(abs.(ax.finallimits[].origin .- last_update_limits[].origin) ./ last_update_limits[].widths)
        change_frac < 0.2 && sleep(Millisecond(500))

        data[]
        last_update_limits[] = lims
        last_update_viewport_widths[] = lon_length, lat_length
    end

    # use colormap if only one layer is supplied
    if dggs isa DGGSArray || (dggs isa DGGSPyramid && args isa Tuple{Symbol})
        filtered_data = filter(x -> !ismissing(x) && !isnan(x), data[])
        cb_limits = (minimum(filtered_data), maximum(filtered_data))

        label = if dggs isa DGGSPyramid && length(args) == 1
            String(args[1])
        else
            DD.label(dggs)
        end

        cb = Colorbar(fig[1, 2], width=10, colormap=:viridis, limits=cb_limits; label=label)
        heatmap!(ax, data, colorrange=cb_limits)
    else
        image!(ax, data)
    end

    return fig
end

end
