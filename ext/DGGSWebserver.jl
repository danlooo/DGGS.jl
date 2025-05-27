module DGGSWebserver

using DGGS
using CairoMakie
using WGLMakie
using Oxygen
using OteraEngine
using HTTP
using DimensionalData

function request_root(collections)
    tmpl = Template("src/html_templates/root.html")
    tmpl(init=Dict(:title => "DGGSExplorer", :collectionIds => keys(collections)))
end

function request_collections_json(collections)
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

function request_collections_html(collections)
    tmpl = Template("src/html_templates/collections.html")
    tmpl(init=Dict(:title => "collections", :collectionIds => keys(collections)))
end

function request_collections(req, collections)
    query_params = queryparams(req)
    f = get(query_params, "f", "json")
    if f == "html"
        return request_collections_html(collections)
    else
        return request_collections_json(collections)
    end
end

function request_collection(req, collectionId, collections)
    query_params = queryparams(req)
    f = get(query_params, "f", "json")
    collection = get(collections, Symbol(collectionId), nothing)

    isnothing(collection) && error("Collection not found: $collectionId")

    if f == "html"
        return request_collection_html(collectionId, collection)
    else
        return request_collection_json(collectionId, collection)
    end
end


function request_collection_html(collectionId, collection::DGGSDataset)
    tmpl = Template("src/html_templates/collection.html")
    tmpl(init=Dict(
        :title => "DGGSExplorer",
        :collectionId => collectionId,
        :layers => keys(layers(collection)),
        :collection => collection,
    ))
end


function request_collection_json(collectionId, collection::DGGSDataset)
    return Dict(
        :id => collectionId
    )
end

function request_collection_map(req, collectionId, collections)
    collection = collections[Symbol(collectionId)]
    fig = CairoMakie.plot(collection, :Red, :Green, :Blue; scale_factor=1 / 255)
    png(fig)
end

function Oxygen.serve(
    collections::Dict{Symbol,T};
    kwargs...
) where {T<:DGGSDataset}
    @get "/" req -> request_root(collections)
    @get "/collections" req -> request_collections(req, collections)
    @get "/collections/{collectionId}" (req, collectionId) -> request_collection(req, collectionId, collections)
    @get "/collections/{collectionId}/map" (req, collectionId) -> request_collection_map(req, collectionId, collections)

    @get "/hello" function ()
        fig = heatmap(rand(50, 50))
        html(fig)
    end

    serve(; kwargs...)
end

end