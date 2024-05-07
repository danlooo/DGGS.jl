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
using ThreadedIterables
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

# must be sorted by dependency
include("types.jl")
include("dggrid.jl")
include("cubes.jl")
include("pyramids.jl")

export
    DGGSArray,
    DGGSPyramid,
    DGGSDatasetPyramid,
    DGGSGridSystem,
    Q2DI,
    transform_points,
    ColorScale,
    query,
    to_dggs_array,
    to_dggs_dataset_pyramid,
    to_geo_cube,
    At,
    BBox
end