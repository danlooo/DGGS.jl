module DGGSMakie

using DGGS
using Makie
using DimensionalData
import DimensionalData as DD
using Makie.GeometryBasics
using Dates

"""
Plot a DGGSArray using Makie.

`resolution_scale` is a value between 0 (calculate no pixel) and 1 (calculate all pixels). Used to reduce plot time.
"""
function Makie.plot(
    dggs_array::DGGSArray,
    args...;
    lon_range=-180:0.25:180,
    lat_range=-90:0.25:90,
    resolution_scale::Real=1,
    kwargs...
)
    length(DGGS.non_spatial_dims(dggs_array)) == 0 || error("DGGSArray must not have any non-spatial dimension")
    0 <= resolution_scale <= 1 || error("resolution_scale must be between 0 and 1")

    fig = Figure()
    ax = Axis(fig[1, 1], limits=(-180, 180, -90, 90), xlabel="longitude [°]", ylabel="latitude [°]")
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
end