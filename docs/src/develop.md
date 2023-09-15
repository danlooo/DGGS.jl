# Developer information


## Create a DGGS

1. Shell call to DGGRID to get grids i.e. `Vector{DgGrid}` at all resolutions. A grid is a ball tree containing cell center points in geographical coordinates
2. Indexing i.e. reorder cells by sorted address, e.g. HIndex
3. Add the data, i.e. create a cell cube based on that grid. Array indices and memory address depends on the calculated index.


Cell center points in geographical space are ordered for fast point queries using `NearestNeighbors.BallTree(;reorder=true)`.
To convert across indices, methods for `Base.getindex` are implemented for index types `<:AbstractIndex`.


```julia
using DGGS
grids = [DgGrid(:isea, 7, :hexagon, resolution) for resolution in 0:4]


indices = [:seqnum]

grid = DgGrid(:isea, 7, :hexagon, 1, :seqnum)
grid[1] # geo coordinates of first cell

grid.data # ball tree of geo coordinates
grid.index 

```
