module DGGS

using DGGRID7_jll
using DimensionalData
using DimensionalData: DimArray, metadata
using Zarr
using NetCDF
using YAXArrays
using DataFrames
using CSV
using ColorSchemes
using FileIO
using IterTools
using ProgressMeter
using Infiltrator
using Statistics
using JSON3
using Makie
using GeometryBasics
using OrderedCollections
using ImageCore
using LinearAlgebra
using FileIO
using MeshIO
using Pkg.Artifacts
using Distributed
using ThreadsX
using Pkg.Artifacts
using FillArrays

# must be sorted by dependency
include("types.jl")
include("dggrid.jl")
include("neighbors.jl")
include("array.jl")
include("layer.jl")
include("pyramid.jl")

export
    DGGSArray,
    DGGSLayer,
    DGGSPyramid,
    Q2DI,
    get_ring,
    getindex,
    transform_points,
    open_dggs_pyramid,
    write_dggs_pyramid,
    to_geo_array,
    to_dggs_array,
    to_dggs_layer,
    to_dggs_pyramid
end