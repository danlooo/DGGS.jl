# The background behind Discrete Global Grid Systems (DGGS)

A Discrete Global Grid Systems (DGGS) tessellate the surface of the earth with hierarchical cells of equal area.
This minimizes distortion and loading time of large geospatial datasets, which is crucial in spatial statistics and building Machine Learning models.

## Why to use a DGGS

- multi resolutions
- chunking and compression
- spherical properties
- equal area

## DGGS Creation

1. Take a platonic solid (e.g. icosahedron)
2. Blow it up so that it s size fits the radius of the earth
3. Chose a rotation of the polyhedron relative to the  (e.g. those used in the [Dymaxion projection](https://en.wikipedia.org/wiki/Dymaxion_map) so that the evrtices of the polyhedrons with high distortions are in the oceans)
4. Tessellate the faces of the polyhedron (e.g. triangles, diamonds, or hexagons). One must introduce 12 pentagons at the vertices of the polyhedron to enable a tesselation of the 3D surface with hexagons.
5. Redo the tessellation with increasing resolutions forming a grid system

![](https://upload.wikimedia.org/wikipedia/commons/thumb/5/53/Dymaxion_projection.png/1920px-Dymaxion_projection.png)
By Justin Kunimune - Own work, Public Domain, https://commons.wikimedia.org/w/index.php?curid=65694588

![](assets/hexagon-children-aperture.png)

![](assets/dggs-distortion.png)