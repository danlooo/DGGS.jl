module DGGSWebserver

using DGGS
using WGLMakie
using Oxygen

function request_root()
    return Dict(
        :links => [
            Dict(
                :rel => "service-doc",
                :href => "/docs",
                :type => "text/html",
                :title => "Swagger documentation",
            ),
        ]
    )
end

function request_collections(collections)
    collections_d = [Dict(
        :id => k,
        :crs => [
            "http://www.opengis.net/def/crs/EPSG/0/4326",
            v.dggsrs
        ],
        :extent => Dict(
            :spatial => Dict(
                :bbox => [DGGS.get_geo_bbox(v) |> x -> [x.X[1], x.Y[1], x.X[2], x.Y[2]]]
            )
        )
    ) for (k, v) in collections]

    return Dict(
        :collections => collections_d
    )
end

function Oxygen.serve(
    collections::Dict{Symbol,T};
    kwargs...
) where {T<:DGGSDataset}
    @get "/" request_root
    @get "/collections" x -> request_collections(collections)

    @get "/hello" function ()
        fig = heatmap(rand(50, 50))
        html(fig)
    end

    serve(; kwargs...)
end

end