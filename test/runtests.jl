using DGGS
using DimensionalData
using YAXArrays
using GLMakie
using Test

base_path = "https://s3.bgc-jena.mpg.de:9000/dggs/sresa1b_ncar_ccsm3-example"

arr = open_array("$base_path/3/tas")
l = open_layer("$base_path/3")
dggs = open_pyramid("$base_path")

@test arr.attrs |> length == 31
@test l.attrs |> length == 19
@test dggs.attrs |> length == 19
@test (setdiff(dggs.attrs, l.attrs) .|> x -> x.first) == ["_DGGS"] # same global attrs expect DGGS level

d = tempname()
write_pyramid(d, dggs)
dggs_2 = open_pyramid(d)
@test dggs.attrs == dggs_2.attrs
@test dggs.levels == dggs_2.levels
@test dggs.bands == dggs_2.bands
rm(d, recursive=true)