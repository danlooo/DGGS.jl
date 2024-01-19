const infobox = document.getElementById("infobox");

/*
 * https://deck.gl/docs/api-reference/core/globe-view
 */
const {
  DeckGL,
  MapView,
  _GlobeView,
  TileLayer,
  BitmapLayer,
  GeoJsonLayer,
  SolidPolygonLayer,
  COORDINATE_SYSTEM,
} = deck;

background_polygon_north = new SolidPolygonLayer({
  id: "background_polygon_north_layer",
  data: [
    [
      [-180, 90],
      [0, 90],
      [180, 90],
      [180, 85],
      [0, 85],
      [-180, 85],
    ],
  ],
  getPolygon: (d) => d,
  stroked: false,
  filled: true,
  getFillColor: [38, 38, 38],
});

background_polygon_south = new SolidPolygonLayer({
  id: "background_polygon_south_layer",
  data: [
    [
      [-180, -85],
      [0, -85],
      [180, -85],
      [180, -90],
      [0, -90],
      [-180, -90],
    ],
  ],
  getPolygon: (d) => d,
  stroked: false,
  filled: true,
  getFillColor: [9, 9, 9],
});

basemap_layer = new TileLayer({
  id: "basemap_layer",
  data: "https://a.basemaps.cartocdn.com/dark_nolabels/{z}/{x}/{y}.png",
  minZoom: 0,
  maxZoom: 19,
  tileSize: 256,

  renderSubLayers: (props) => {
    const {
      bbox: { west, south, east, north },
    } = props.tile;

    return new BitmapLayer(props, {
      data: null,
      image: props.data,
      _imageCoordinateSystem: COORDINATE_SYSTEM.CARTESIAN,
      bounds: [west, south, east, north],
    });
  },
});

data_layer = new TileLayer({
  id: "data_layer",
  data: "collections/data%252Fmodis-ndvi.dggs/tiles/{z}/{x}/{y}",
  minZoom: 0,
  maxZoom: 19,
  tileSize: 256,
  // only one scheduling queue to allow aborting not required tiles after pan and zoom
  // see https://github.com/visgl/deck.gl/issues/4429
  maxRequests: 1,
  maxOngoingRequests: 4,

  renderSubLayers: (props) => {
    const {
      bbox: { west, south, east, north },
    } = props.tile;

    return new BitmapLayer(props, {
      data: null,
      image: props.data,
      _imageCoordinateSystem: COORDINATE_SYSTEM.CARTESIAN,
      bounds: [west, south, east, north],
    });
  },

  onViewStateChange: ({ viewState }) => {
    infobox.innerHTML = viewState.zoom;
  },
});

coastline_layer = new GeoJsonLayer({
  id: "coastline_layer",
  data: "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_coastline.geojson",
  getLineColor: [255, 255, 255],
  lineWidthMinPixels: 1,
  lineWidthMaxPixels: 1,
});

new DeckGL({
  // views: new _GlobeView({
  //   resolution: 10,
  // }),
  views: new MapView(),
  initialViewState: {
    longitude: 0,
    latitude: 0,
    zoom: 1,
    minZoom: 0,
    maxZoom: 20,
  },
  controller: true,
  layers: [
    background_polygon_north,
    background_polygon_south,
    basemap_layer,
    data_layer,
    coastline_layer,
  ],
});
