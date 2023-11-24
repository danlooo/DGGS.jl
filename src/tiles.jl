"Rectangular bounding box in geographical space"
struct BBox{T<:Real}
    lon_min::T
    lon_max::T
    lat_min::T
    lat_max::T
end

"Geographical bounding box given a xyz tile"
function BBox(x, y, z)
    #@see https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames

    n = 2^z
    lon_min = x / n * 360.0 - 180.0
    lat_min = atan(sinh(π * (1 - 2 * (y + 1) / n))) |> rad2deg

    lon_max = (x + 1) / n * 360.0 - 180.0
    lat_max = atan(sinh(π * (1 - 2 * y / n))) |> rad2deg

    return BBox(lon_min, lon_max, lat_min, lat_max)
end

struct ColorScale{T<:Real}
    schema::ColorScheme
    min_value::T
    max_value::T
end


function color(value::Real, color_scheme::ColorScheme; null_color=RGBA{Float64}(0, 0, 0, 0))
    isnan(value) && return null_color
    ismissing(value) && return null_color
    return color_scale.schema[value] |> RGBA
end

