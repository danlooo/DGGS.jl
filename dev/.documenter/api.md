<details class='jldocstring custom-block' open>
<summary><a id='Base.getindex-Tuple{DGGSArray, Integer, Integer, Integer, Vararg{Any}}' href='#Base.getindex-Tuple{DGGSArray, Integer, Integer, Integer, Vararg{Any}}'><span class="jlbinding">Base.getindex</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



get a cell of a DGGSArray


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L59" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.getindex-Tuple{DGGSArray, Q2DI, Integer}' href='#Base.getindex-Tuple{DGGSArray, Q2DI, Integer}'><span class="jlbinding">Base.getindex</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



get a ring of a DGGArray


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L65" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.getindex-Tuple{DGGSArray, Q2DI}' href='#Base.getindex-Tuple{DGGSArray, Q2DI}'><span class="jlbinding">Base.getindex</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



get a cell of a DGGSArray


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L55" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.getindex-Tuple{DGGSArray}' href='#Base.getindex-Tuple{DGGSArray}'><span class="jlbinding">Base.getindex</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



filter any dimension of a DGGSArray


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.getindex-Union{Tuple{R}, Tuple{DGGSArray, Q2DI, UnitRange{R}}} where R<:Integer' href='#Base.getindex-Union{Tuple{R}, Tuple{DGGSArray, Q2DI, UnitRange{R}}} where R<:Integer'><span class="jlbinding">Base.getindex</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



get a disk of a DGGArray


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L78" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS._transform_points-Union{Tuple{V}, Tuple{U}, Tuple{AbstractArray{Tuple{U, V}, 1}, Any}} where {U<:Real, V<:Real}' href='#DGGS._transform_points-Union{Tuple{V}, Tuple{U}, Tuple{AbstractArray{Tuple{U, V}, 1}, Any}} where {U<:Real, V<:Real}'><span class="jlbinding">DGGS._transform_points</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Transforms Vector of (lon,lat) coords to DGGRID indices


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/dggrid.jl#L55" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.call_dggrid-Tuple{Dict}' href='#DGGS.call_dggrid-Tuple{Dict}'><span class="jlbinding">DGGS.call_dggrid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Execute sytem call of DGGRID binary


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/dggrid.jl#L4-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.clip-Tuple{UnitRange, Any}' href='#DGGS.clip-Tuple{UnitRange, Any}'><span class="jlbinding">DGGS.clip</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



clip range that may lie in padding to quad boundaries


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/neighbors.jl#L29-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.disk_mask-Tuple{Integer}' href='#DGGS.disk_mask-Tuple{Integer}'><span class="jlbinding">DGGS.disk_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Create a mask of cells arround a center cell within a radius `disk_size`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/neighbors.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.filter_null-Tuple{Any}' href='#DGGS.filter_null-Tuple{Any}'><span class="jlbinding">DGGS.filter_null</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Apply function f after filtering of missing and NAN values


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L167" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.plot_native-Tuple{DGGSArray}' href='#DGGS.plot_native-Tuple{DGGSArray}'><span class="jlbinding">DGGS.plot_native</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Plot a DGGSArray nativeley on a icosahedron


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L517-L519" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.ring_mask-Tuple{Integer}' href='#DGGS.ring_mask-Tuple{Integer}'><span class="jlbinding">DGGS.ring_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Create a mask of cells having the same distance `ring_size` to a center cell


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/neighbors.jl#L10-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.to_dggs_array-Tuple{DimensionalData.AbstractDimArray, Integer}' href='#DGGS.to_dggs_array-Tuple{DimensionalData.AbstractDimArray, Integer}'><span class="jlbinding">DGGS.to_dggs_array</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Convert a geographic lat/lon raster to a DGGAArray

agg_type: `:convert` will return a Float64 array and `:round` will keep the element type that might loose precision


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/array.jl#L180-L184" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.to_dggs_layer-Tuple{YAXArrays.Datasets.Dataset, Integer}' href='#DGGS.to_dggs_layer-Tuple{YAXArrays.Datasets.Dataset, Integer}'><span class="jlbinding">DGGS.to_dggs_layer</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Transforms a `YAXArrays.Dataset` in geographic lat/lon ratser to a DGGSLayer at agiven layer


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/layer.jl#L97-L99" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.transform_points-Union{Tuple{B}, Tuple{A}, Tuple{AbstractVector{A}, AbstractVector{B}, Integer}} where {A<:Real, B<:Real}' href='#DGGS.transform_points-Union{Tuple{B}, Tuple{A}, Tuple{AbstractVector{A}, AbstractVector{B}, Integer}} where {A<:Real, B<:Real}'><span class="jlbinding">DGGS.transform_points</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



chunk_size_points: number of points (e.g. pixels) to transform in one block (task of a thread)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/dggrid.jl#L144-L146" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='DGGS.width-Tuple{Integer}' href='#DGGS.width-Tuple{Integer}'><span class="jlbinding">DGGS.width</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



position of last row or column in a quad matrix of that level


<Badge type="info" class="source-link" text="source"><a href="https://github.com/danlooo/DGGS.jl/blob/8a9d996cc4830be9d91d4ee967ff86d960ca5d05/src/pyramid.jl#L93" target="_blank" rel="noreferrer">source</a></Badge>

</details>

