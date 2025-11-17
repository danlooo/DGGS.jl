module DGGS

import Proj
using Base.Threads
using YAXArrays
using DimensionalData
import DimensionalData as DD
using Statistics
using CoordinateTransformations
using GeometryBasics
using Dates
using Printf
using Infiltrator
using Extents
using OrderedCollections
using DiskArrayTools

include("types.jl")
include("cells.jl")
include("arrays.jl")
include("datasets.jl")
include("pyramids.jl")

const transformations = Channel{Proj.Transformation}(Inf)
const inv_transformations = Channel{Proj.Transformation}(Inf)
const threads_ready = Ref(false)

export Cell, DGGSArray, DGGSDatase, Cell, DGGSArray, DGGSDataset, DGGSPyramid
export to_cell, to_geo, to_dggs_pyramid, to_dggs_dataset, to_dggs_array, to_geo_dataset, to_geo_array
export open_dggs_array, open_dggs_dataset, open_dggs_pyramid
export save_dggs_array, save_dggs_dataset, save_dggs_pyramid
end