function run_webserver(; kwargs...)
  @swagger """
  /collections/{path}/tiles/{z}/{x}/{y}:
    get:
      description: Calculate a XYZ tile. DGGS datacube my be filtered using URL query parameters
      parameters:
        - name: path
          in: path
          required: true
          description: URL encoded path of the DGGS data cube
          schema:
            type : string
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
        - name: min_val
          in: query
          required: false
          description: Lowest value used as the lower bound for the color scale.
          schema:
            type : number
        - name: max_val
          in: query
          required: false
          description: Highest value used as the upper bound for the color scale.
          schema:
            type : number
      responses:
        '200':
          description: Successfully returned the tile
  """
  @get "/collections/{path}/tiles/{z}/{x}/{y}" function (req::HTTP.Request, path::String, z::Int, x::Int, y::Int)
    dggs = path |> HTTP.unescapeuri |> HTTP.unescapeuri |> GridSystem

    params = HTTP.queryparams(req)
    min_val = get(params, "min_val", "0") |> x -> parse(Float64, x)
    max_val = get(params, "max_val", "1") |> x -> parse(Float64, x)

    query_d = filter(((k, v),) -> k ∉ ["min_val", "max_val"], params)
    query_str = ""
    for (k, v) in query_d
      query_str *= "$k=$v"
    end

    color_scale = ColorScale(ColorSchemes.viridis, min_val, max_val)
    tile = calculate_tile(dggs, color_scale, x, y, z; query_str=query_str, cache_path="data/cache_xyz_to_q2di")
    response_headers = [
      "Content-Type" => "image/png",
      # "cache-control" => "max-age=23117, stale-while-revalidate=604800, stale-if-error=604800"
    ]
    response = HTTP.Response(200, response_headers, tile)
    return response
  end

  @get "/collections/{path}" function (req::HTTP.Request, path::String)
    path = HTTP.unescapeuri(path)
    JSON3.read("$(path)/.zattrs")
  end


  # TODO: Use Artifacts, see https://github.com/JuliaPackaging/ArtifactUtils.jl to upload directory to github
  dynamicfiles("/home/dloos/prj/DGGS.jl/src/assets/www", "/")

  info = Dict("title" => "DGGSexplorer API", "version" => "1.0.0")
  openApi = OpenAPI("3.0", info)
  swagger_document = build(openApi)
  mergeschema(swagger_document)

  serveparallel(; kwargs...)
end