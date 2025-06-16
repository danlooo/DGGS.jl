```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: DGGS.jl
  text: Discrete Global Grid System for Julia
  tagline: Less distorted geospatial data cubes
  image:
    src: icon.drawio.svg
    alt: VitePress
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/danlooo/DGGS.jl
    - theme: alt
      text: API reference
      link: /api
features:
  - title: Preview
    details: Explore the preliminary version of <a href="https://github.com/danlooo/DGGS.jl/tree/pentacube">DGGS.jl 2.0</a> using the faster pentacube index build with PROJ. We are also developing <a href="http://dggs.fairsendd.eodchosting.eu">DGGSexplorer</a>: A web server to view DGGS data cubes in a browser. 
  - title: Low distortions
    details: This package is based on <a href="https://github.com/sahrk/DGGRID">DGGRID</a> which has the lowest distortion in shape and area compared to other DGGS. 
  - title: DGGS native data cubes
    details: The geospatial data is directly stored along dimensions of a DGGS index to optimize space and time needed to store, load, and process the data.
  - title: Features
    details: Conversion of coordinates and data between geographical and DGGS spaces, data visualization, Neighbor retrieval
```
