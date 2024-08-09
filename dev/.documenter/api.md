<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Base.getindex-Tuple{DGGSArray, Q2DI, Integer}' href='#Base.getindex-Tuple{DGGSArray, Q2DI, Integer}'>#</a>&nbsp;<b><u>Base.getindex</u></b> &mdash; <i>Method</i>.




get a ring of a DGGArray


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/array.jl#L65)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Base.getindex-Tuple{DGGSArray, Q2DI}' href='#Base.getindex-Tuple{DGGSArray, Q2DI}'>#</a>&nbsp;<b><u>Base.getindex</u></b> &mdash; <i>Method</i>.




get a cell of a DGGSArray


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/array.jl#L55)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Base.getindex-Tuple{DGGSArray}' href='#Base.getindex-Tuple{DGGSArray}'>#</a>&nbsp;<b><u>Base.getindex</u></b> &mdash; <i>Method</i>.




filter any dimension of a DGGSArray


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/array.jl#L31)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Base.getindex-Union{Tuple{R}, Tuple{DGGSArray, Q2DI, UnitRange{R}}} where R<:Integer' href='#Base.getindex-Union{Tuple{R}, Tuple{DGGSArray, Q2DI, UnitRange{R}}} where R<:Integer'>#</a>&nbsp;<b><u>Base.getindex</u></b> &mdash; <i>Method</i>.




get a disk of a DGGArray


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/array.jl#L78)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Base.getindex-Union{Tuple{V}, Tuple{U}, Tuple{T}, Tuple{DGGSArray, T, U, V}} where {T<:Integer, U<:Integer, V<:Integer}' href='#Base.getindex-Union{Tuple{V}, Tuple{U}, Tuple{T}, Tuple{DGGSArray, T, U, V}} where {T<:Integer, U<:Integer, V<:Integer}'>#</a>&nbsp;<b><u>Base.getindex</u></b> &mdash; <i>Method</i>.




get a cell of a DGGSArray


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/array.jl#L59)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS._transform_points-Union{Tuple{V}, Tuple{U}, Tuple{AbstractArray{Tuple{U, V}, 1}, Any}} where {U<:Real, V<:Real}' href='#DGGS._transform_points-Union{Tuple{V}, Tuple{U}, Tuple{AbstractArray{Tuple{U, V}, 1}, Any}} where {U<:Real, V<:Real}'>#</a>&nbsp;<b><u>DGGS._transform_points</u></b> &mdash; <i>Method</i>.




Transforms Vector of (lon,lat) coords to DGGRID indices


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/dggrid.jl#L55)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.aggregate_dggs_array-Tuple{Any, Any, DGGSArray}' href='#DGGS.aggregate_dggs_array-Tuple{Any, Any, DGGSArray}'>#</a>&nbsp;<b><u>DGGS.aggregate_dggs_array</u></b> &mdash; <i>Method</i>.




Spatial hexagonal convolution in Q2DI index space matching levels of DGGRID ISEA4H grids


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/pyramid.jl#L216-L218)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.call_dggrid-Tuple{Dict}' href='#DGGS.call_dggrid-Tuple{Dict}'>#</a>&nbsp;<b><u>DGGS.call_dggrid</u></b> &mdash; <i>Method</i>.




Execute sytem call of DGGRID binary


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/dggrid.jl#L4-L6)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.clip-Tuple{UnitRange, Any}' href='#DGGS.clip-Tuple{UnitRange, Any}'>#</a>&nbsp;<b><u>DGGS.clip</u></b> &mdash; <i>Method</i>.




clip range that may lie in padding to quad boundaries


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/neighbors.jl#L29-L31)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.disk_mask-Tuple{Integer}' href='#DGGS.disk_mask-Tuple{Integer}'>#</a>&nbsp;<b><u>DGGS.disk_mask</u></b> &mdash; <i>Method</i>.




Create a mask of cells arround a center cell within a radius `disk_size`


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/neighbors.jl#L1-L3)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.filter_null-Tuple{Any}' href='#DGGS.filter_null-Tuple{Any}'>#</a>&nbsp;<b><u>DGGS.filter_null</u></b> &mdash; <i>Method</i>.




Apply function f after filtering of missing and NAN values


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/array.jl#L167)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.plot_native-Tuple{DGGSArray}' href='#DGGS.plot_native-Tuple{DGGSArray}'>#</a>&nbsp;<b><u>DGGS.plot_native</u></b> &mdash; <i>Method</i>.




Plot a DGGSArray nativeley on a icosahedron


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/array.jl#L506-L508)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.ring_mask-Tuple{Integer}' href='#DGGS.ring_mask-Tuple{Integer}'>#</a>&nbsp;<b><u>DGGS.ring_mask</u></b> &mdash; <i>Method</i>.




Create a mask of cells having the same distance `ring_size` to a center cell


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/neighbors.jl#L10-L12)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.to_dggs_layer-Tuple{YAXArrays.Datasets.Dataset, Integer}' href='#DGGS.to_dggs_layer-Tuple{YAXArrays.Datasets.Dataset, Integer}'>#</a>&nbsp;<b><u>DGGS.to_dggs_layer</u></b> &mdash; <i>Method</i>.




Transforms a `YAXArrays.Dataset` in geographic lat/lon ratser to a DGGSLayer at agiven layer


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/layer.jl#L99-L101)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.transform_points-Union{Tuple{B}, Tuple{A}, Tuple{AbstractVector{A}, AbstractVector{B}, Integer}} where {A<:Real, B<:Real}' href='#DGGS.transform_points-Union{Tuple{B}, Tuple{A}, Tuple{AbstractVector{A}, AbstractVector{B}, Integer}} where {A<:Real, B<:Real}'>#</a>&nbsp;<b><u>DGGS.transform_points</u></b> &mdash; <i>Method</i>.




chunk_size_points: number of points (e.g. pixels) to transform in one block (task of a thread)


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/dggrid.jl#L139-L141)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='DGGS.width-Tuple{Integer}' href='#DGGS.width-Tuple{Integer}'>#</a>&nbsp;<b><u>DGGS.width</u></b> &mdash; <i>Method</i>.




position of last row or column in a quad matrix of that level


[source](https://github.com/danlooo/DGGS.jl/blob/a21aa320207cfd8b52ce7b50eeb537a8afba76a5/src/pyramid.jl#L100)

</div>
<br>
