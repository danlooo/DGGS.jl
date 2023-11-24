/* **** Leaflet **** */

var base_layer = L.tileLayer(
  "https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png",
  {
    minZoom: 0,
    maxZoom: 16,
  }
);

var data_layer = L.tileLayer("tile/{z}/{x}/{y}/tile.png", {
  layer: "simple",
  minZoom: 0,
  maxZoom: 16,
});

var map = L.map("map", {
  center: [0, 0],
  zoom: 3,
  minZoom: 0,
  maxZoom: 7,
  layers: [base_layer, data_layer],
});

fetch(
  "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_coastline.geojson"
)
  .then((x) => x.json())
  .then((x) =>
    L.geoJSON(x, { style: { color: "#ffffff", weight: 0.7 } }).addTo(map)
  );

var layerControl = L.control.layers([base_layer, data_layer], {}).addTo(map);
