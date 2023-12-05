/*
 * https://deck.gl/docs/api-reference/core/globe-view
 */
const {
  DeckGL,
  _GlobeView,
  TileLayer,
  BitmapLayer,
  GeoJsonLayer,
  COORDINATE_SYSTEM,
} = deck;

new DeckGL({
  views: new _GlobeView({
    resolution: 10,
  }),
  initialViewState: {
    longitude: 0,
    latitude: 0,
    zoom: 1,
    minZoom: 0,
    maxZoom: 20,
  },
  controller: true,

  layers: [
    new TileLayer({
      id: "basemap_layer",
      data: "https://a.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png",
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
          bounds: [west, south, east, north],
        });
      },
    }),
    new TileLayer({
      id: "data_layer",
      data: "tile/{z}/{x}/{y}/tile.png",
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
    }),
    new GeoJsonLayer({
      id: "coastline_layer",
      data: "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_coastline.geojson",
      getLineWidth: 10,
    }),
  ],
});
