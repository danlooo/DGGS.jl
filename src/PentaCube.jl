module PentaCube

using Infiltrator
import Proj
using Base.Threads
using YAXArrays
using DimensionalData
using Statistics

include("types.jl")
include("cells.jl")
include("arrays.jl")

function __init__()
    @info "PentaCube initialized with $(nthreads()) threads"
    global transformations = [Proj.Transformation(crs_geo, crs_isea; ctx=Proj.proj_context_create()) for _ in 1:nthreads()]
    global inv_transformations = [Proj.Transformation(crs_isea, crs_geo; ctx=Proj.proj_context_create()) for _ in 1:nthreads()]
end

export to_cell, to_geo, to_dggs_array, Cell
end