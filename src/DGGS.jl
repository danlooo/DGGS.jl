module DGGS

import Proj
using Base.Threads
using YAXArrays
using DimensionalData
import DimensionalData as DD
using Statistics
using CoordinateTransformations

include("types.jl")
include("cells.jl")
include("arrays.jl")

const transformations = Channel{Proj.Transformation}(Inf)
const inv_transformations = Channel{Proj.Transformation}(Inf)
const threads_ready = Ref(false)

function __init__()
    for _ in 1:Threads.nthreads()
        put!(transformations, Proj.Transformation(crs_geo, crs_isea; ctx=Proj.proj_context_create()))
        put!(inv_transformations, Proj.Transformation(crs_isea, crs_geo; ctx=Proj.proj_context_create()))
    end
    @info "DGGS.jl initialized with $(nthreads()) threads"
end

export to_cell, to_geo, to_dggs_array, to_geo_array, Cell
end