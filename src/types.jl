struct ColorScale{T<:Real}
    schema::ColorScheme
    min_value::T
    max_value::T
end

struct Q2DI
    n
    i
    j
end

struct CellCube
    data::YAXArray
    level::Int8
end

struct GeoCube
    data::YAXArray

    function GeoCube(data)
        :lon in propertynames(data) || error("Axis with name :lon must be present")
        :lat in propertynames(data) || error("Axis with name :lat must be present")
        -180 <= minimum(data.lon) <= maximum(data.lon) <= 180 || error("All longitudes must be within [-180, 180]")
        -90 <= minimum(data.lat) <= maximum(data.lat) <= 90 || error("All latitudes must be within [-90, 90]")

        new(data)
    end
end

struct GridSystem
    data::Dict{Int,CellCube} # some levels might be skipped

    function GridSystem(data)
        map(x -> x.data.axes .|> name, values(data)) |> allequal || error("Same dimensions must be used at all levels.")
        new(data)
    end
end