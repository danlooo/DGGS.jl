module DGGSMakie

using DGGS
using Makie
using DimensionalData
using Makie.GeometryBasics

function Makie.plot(dggs_array::DGGSArray, args...; kwargs...)
    if length(DGGS.non_spatial_dims(dggs_array)) == 0
        plot_single_var(dggs_array, args...; kwargs...)
    else
        plot_multi_var(dggs_array, args...; kwargs...)
    end
end

function plot_single_var(
    dggs_array::DGGSArray, args...;
    lon_range=-180:0.25:180, lat_range=-90:0.25:90,
    resolution=1000,
    kwargs...
)
    length(DGGS.non_spatial_dims(dggs_array)) == 0 || error("DGGSArray must not have any non-spatial dimension")

    fig = Figure()
    ax = Axis(fig[1, 1], limits=(-180, 180, -90, 90))

    data = Observable(to_geo_array(dggs_array, lon_range, lat_range))

    # update plot after drag
    register_interaction!(ax, :my_mouse_interaction) do event::MouseEvent, axis
        if event.type == Makie.MouseEventTypes.rightdragstop
            lims = axis.finallimits[]
            lat_resolution = Int(round(lims.widths[2] / lims.widths[1] * resolution))
            lon_range = range(lims.origin[1], lims.origin[1] + lims.widths[1]; length=resolution)
            lat_range = range(lims.origin[2], lims.origin[2] + lims.widths[2]; length=lat_resolution)

            data[] = to_geo_array(dggs_array, lon_range, lat_range)
        end
    end

    # TODO: add update plot after zoom  

    # enforce zoom to be inside of lat/lon limits
    on(ax.finallimits) do fl
        if abs(fl.origin[1]) + abs(fl.widths[1]) > 180 || abs(fl.origin[2]) + abs(fl.widths[2]) > 90
            ax.targetlimits[] = HyperRectangle{2,Float32}([-180, -90], [360, 180])
        end
    end

    heatmap!(ax, data)
    fig
end

function plot_multi_var(dggs_array::DGGSArray, args...; kwargs...)
    #TODO: implement
    fig = Figure()
    ax = Axis(fig[1, 1:2])
    s1 = Slider(fig[2, 1], range=0.1:0.1:10, startvalue=3)
    data = lift(s1.value) do v
        1:v
    end
    p = scatter!(fig[1, 1:2], data, markersize=s1.value)

    fig
end

end