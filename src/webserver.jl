module WebServer
using Oxygen

function run()
  @swagger """
  /tile/{z}/{x}/{y}/tile.png:
    get:
      description: Calculate a XYZ tile
      parameters:
        - name: x
          in: path
          required: true
          description: Column [OSM](https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames)
          schema:
            type : number
        - name: y
          in: path
          required: true
          description: Row see [OSM](https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames)
          schema:
            type : number
        - name: z
          in: path
          required: true
          description: Zoom level see [OSM](https://wiki.openstreetmap.org/wiki/Zoom_levels)
          schema:
            type : number
      responses:
        '200':
          description: Successfully returned an number.
  """
  @get "/tile/{z}/{x}/{y}/tile.png" function (req::HTTP.Request, z::Int, x::Int, y::Int)
    params = HTTP.queryparams(req)

    tile = calculate_tile(dggs, cell_ids, color_scale, x, y, z)
    response_headers = [
      "Content-Type" => "image/png",
      # "cache-control" => "max-age=23117, stale-while-revalidate=604800, stale-if-error=604800"
    ]
    response = HTTP.Response(200, response_headers, tile)
    return response
  end
  dynamicfiles("www", "/")

  info = Dict("title" => "DGGSexplorer API", "version" => "1.0.0")
  openApi = OpenAPI("3.0", info)
  swagger_document = build(openApi)
  mergeschema(swagger_document)

  serve(; port=8089)
end

end