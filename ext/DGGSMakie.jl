module DGGSMakie

using DGGS
using Makie
using DimensionalData
import DimensionalData as DD
using Makie.GeometryBasics
using Dates


using Infiltrator

"""
Plot a DGGSArray using Makie.

`resolution_scale` is a value between 0 (calculate no pixel) and 1 (calculate all pixels). Used to reduce plot time.
"""
function Makie.plot(
    dggs_array::DGGSArray,
    args...;
    extent=dggs_array.bbox,
    resolution_scale::Real=1,
    kwargs...
)
    length(DGGS.non_spatial_dims(dggs_array)) == 0 || error("DGGSArray must not have any non-spatial dimension")
    0 <= resolution_scale <= 1 || error("resolution_scale must be between 0 and 1")

    # start zoom on full extent
    lon_range = range(extent.X..., 100)
    lat_range = range(extent.Y..., 100)

    fig = Figure()
    ax = Axis(fig[1, 1], limits=(extent.X..., extent.Y...), xlabel="longitude [°]", ylabel="latitude [°]")
    ax.aspect = DataAspect()

    data = Observable(to_geo_array(dggs_array, lon_range, lat_range))

    last_update_limits = Observable(ax.finallimits[])
    last_update_viewport_widths = Observable(fig.scene.viewport[].widths)

    filtered_data = filter(x -> !ismissing(x) && !isnan(x), data[])
    cb_limits = (minimum(filtered_data), maximum(filtered_data))
    cb = Colorbar(fig[1, 2], width=10, colormap=:viridis, limits=cb_limits; label=DD.label(dggs_array))

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

        # get limits observable once. It might change during the update
        lims = ax.finallimits[]
        lon_min, lat_min = lims.origin
        lon_max, lat_max = lims.origin .+ lims.widths

        lon_min = clamp(lon_min, -180, 180)
        lon_max = clamp(lon_max, -180, 180)
        lat_min = clamp(lat_min, -90, 90)
        lat_max = clamp(lat_max, -90, 90)

        lon_length, lat_length = fig.scene.viewport[].widths
        lon_range = range(lon_min, lon_max, length=lon_length * resolution_scale)
        lat_range = range(lat_min, lat_max, length=lat_length * resolution_scale)
        data[] = to_geo_array(dggs_array, lon_range, lat_range)

        last_update_limits[] = lims
        last_update_viewport_widths[] = lon_length, lat_length
    end

    heatmap!(ax, data, colorrange=cb_limits)
    fig
end

""
function Makie.plot(
    ds::DGGSDataset,
    red_layer::Symbol,
    green_layer::Symbol,
    blue_layer::Symbol,
    args...;
    extent=ds.bbox,
    resolution_scale::Real=1,
    scale_factor::Real=1,
    offset::Real=0,
    kwargs...
)
    red_layer in keys(ds) || error("Red layer $(red_layer) not found in dataset")
    green_layer in keys(ds) || error("Green layer $(green_layer) not found in dataset")
    blue_layer in keys(ds) || error("Blue layer $(blue_layer) not found in dataset")

    # start with zoom to full extent
    lon_dim = Observable(X(range(extent.X..., 500)))
    lat_dim = Observable(Y(range(extent.Y..., 500)))

    data = @lift begin
        ds_rgb = DGGSDataset(getproperty(ds, red_layer), getproperty(ds, green_layer), getproperty(ds, blue_layer))
        geo_ds = to_geo_dataset(ds_rgb, $lon_dim, $lat_dim)

        picture = map(CartesianIndices(($lon_dim, $lat_dim))) do i
            r = getproperty(geo_ds, red_layer)[i] * scale_factor + offset
            g = getproperty(geo_ds, green_layer)[i] * scale_factor + offset
            b = getproperty(geo_ds, blue_layer)[i] * scale_factor + offset

            if ismissing(r) || ismissing(g) || ismissing(b)
                Makie.RGBf(1, 1, 1)
            else
                Makie.RGBf(r, g, b)
            end
        end
        picture
    end

    fig = Figure()
    ax = Axis(fig[1, 1], limits=(extent.X..., extent.Y...), xlabel="longitude [°]", ylabel="latitude [°]")
    ax.aspect = DataAspect()

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

        # get limits observable once. It might change during the update
        lims = ax.finallimits[]
        lon_min, lat_min = lims.origin
        lon_max, lat_max = lims.origin .+ lims.widths

        lon_min = clamp(lon_min, -180, 180)
        lon_max = clamp(lon_max, -180, 180)
        lat_min = clamp(lat_min, -90, 90)
        lat_max = clamp(lat_max, -90, 90)

        lon_length, lat_length = fig.scene.viewport[].widths
        lon_range = range(lon_min, lon_max, length=lon_length * resolution_scale)
        lat_range = range(lat_min, lat_max, length=lat_length * resolution_scale)
        lon_dim[] = X(lon_range)
        lat_dim[] = Y(lat_range)
        data[] # hold loop during plotting the new frame

        last_update_limits[] = lims
        last_update_viewport_widths[] = lon_length, lat_length
    end

    image!(ax, data)
    fig
end
end
