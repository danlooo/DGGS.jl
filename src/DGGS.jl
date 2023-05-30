module DGGS

include("grid.jl")
include("dggrid.jl")

export GridSpec, Grid, toyGrid, PresetGridSpecs, cell_centers, cell_boundaries, dg_call, cell_name, geo_coords
end
