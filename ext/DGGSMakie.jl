module DGGSMakie

using DGGS
using Makie
using DimensionalData
import DimensionalData as DD
using Makie.GeometryBasics
using Dates

function Makie.plot(
    dggs_array::DGGSArray, args...;
    lon_range=-180:0.25:180,
    lat_range=-90:0.25:90,
    resolution=1000,
    update_interval=Millisecond(1500),
    kwargs...
)
    length(DGGS.non_spatial_dims(dggs_array)) == 0 || error("DGGSArray must not have any non-spatial dimension")

    fig = Figure()
    axis = Axis(fig[1, 1], limits=(-180, 180, -90, 90))

    data = Observable(to_geo_array(dggs_array, lon_range, lat_range))
    last_update_limits = Observable(axis.finallimits[])

    filtered_data = filter(x -> !ismissing(x) && !isnan(x), data[])
    cb_limits = (minimum(filtered_data), maximum(filtered_data))
    cb = Colorbar(fig[1, 2], width=10, colormap=:viridis, limits=cb_limits; label=DD.label(dggs_array))

    # update plot after every update interval
    update_time = Observable(now())
    @async while true
        sleep(update_interval)
        update_time[] = now()
    end

    on(update_time) do current_time
        # skip update if limits did not change
        axis.finallimits[] == last_update_limits[] && return

        lon_min, lat_min = axis.finallimits[].origin
        lon_max, lat_max = axis.finallimits[].origin .+ axis.finallimits[].widths

        lims = axis.finallimits[]
        lat_resolution = Int(round(lims.widths[2] / lims.widths[1] * resolution))
        lon_range = range(lon_min, lon_max, length=resolution)
        lat_range = range(lat_min, lat_max, length=lat_resolution)
        data[] = to_geo_array(dggs_array, lon_range, lat_range)
        last_update_limits[] = axis.finallimits[]
    end

    on(axis.finallimits) do lims
        # enforce zoom to be inside of lat/lon limits
        if abs(lims.origin[1]) + abs(lims.widths[1]) > 180 || abs(lims.origin[2]) + abs(lims.widths[2]) > 90
            axis.targetlimits[] = HyperRectangle{2,Float32}([-180, -90], [360, 180])
            axis.finallimits[] = axis.targetlimits[]
        end
    end

    heatmap!(axis, data, colorrange=cb_limits)
    fig
end
end