import{_ as a,c as n,a5 as i,o as p}from"./chunks/framework.ChCwx2Gq.js";const t="/DGGS.jl/dev/assets/rings-disks.DQjhkQvs.png",r=JSON.parse('{"title":"Select","description":"","frontmatter":{},"headers":[],"relativePath":"select.md","filePath":"select.md","lastUpdated":null}'),l={name:"select.md"};function e(o,s,h,d,g,c){return p(),n("div",null,s[0]||(s[0]=[i(`<h1 id="select" tabindex="-1">Select <a class="header-anchor" href="#select" aria-label="Permalink to &quot;Select&quot;">​</a></h1><h2 id="Select-arrays" tabindex="-1">Select arrays <a class="header-anchor" href="#Select-arrays" aria-label="Permalink to &quot;Select arrays {#Select-arrays}&quot;">​</a></h2><p>Open a array for testing:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DGGS</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> open_dggs_pyramid</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;https://s3.bgc-jena.mpg.de:9000/dggs/datasets/modis&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>DGGSPyramid</span></span>
<span class="line"><span>DGGS: DGGRID ISEA4H Q2DI ⬢</span></span>
<span class="line"><span>Levels: [2, 3, 4, 5, 6, 7, 8, 9, 10]</span></span>
<span class="line"><span>Non spatial axes:</span></span>
<span class="line"><span>  Ti 216 Dates.DateTime points</span></span>
<span class="line"><span>Arrays:</span></span>
<span class="line"><span>  lst (:Ti) K Union{Missing, Float32} </span></span>
<span class="line"><span>  ndvi (:Ti) NDVI Union{Missing, Float32}</span></span></code></pre></div><p>Select a ndvi at a given spatial resolution level:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>DGGSLayer{5}</span></span>
<span class="line"><span>DGGS:	DGGRID ISEA4H Q2DI ⬢ at level 5</span></span>
<span class="line"><span>Non spatial axes:</span></span>
<span class="line"><span>  Ti 216 Dates.DateTime points</span></span>
<span class="line"><span>Arrays:</span></span>
<span class="line"><span>  lst (:Ti) K Union{Missing, Float32} </span></span>
<span class="line"><span>  ndvi (:Ti) NDVI Union{Missing, Float32}</span></span></code></pre></div><p>Select an array by its id:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ndvi</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>DGGSArray{Union{Missing, Float32}, 5}</span></span>
<span class="line"><span>Name:		ndvi</span></span>
<span class="line"><span>Units:		NDVI</span></span>
<span class="line"><span>DGGS:		DGGRID ISEA4H Q2DI ⬢ at level 5</span></span>
<span class="line"><span>Attributes:	18</span></span>
<span class="line"><span>Non spatial axes:</span></span>
<span class="line"><span>  Ti 216 Dates.DateTime points</span></span></code></pre></div><p>Additional filtering by any non-spatial axes e.g. <code>Time</code> still results in a <code>DGGSArray</code>:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ndvi[Time</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>DGGSArray{Union{Missing, Float32}, 5}</span></span>
<span class="line"><span>Name:		NDVI</span></span>
<span class="line"><span>Units:		NDVI</span></span>
<span class="line"><span>DGGS:		DGGRID ISEA4H Q2DI ⬢ at level 5</span></span>
<span class="line"><span>Attributes:	18</span></span>
<span class="line"><span>Non spatial axes:</span></span></code></pre></div><p>Further filtering will return a <code>YAXArray</code> instead:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ndvi[Time</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][q2di_n </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭───────────────────────────────────────────╮</span></span>
<span class="line"><span>│ 16×16 YAXArray{Union{Missing, Float32},2} │</span></span>
<span class="line"><span>├───────────────────────────────────────────┴───────────── dims ┐</span></span>
<span class="line"><span>  ↓ q2di_i Sampled{Int64} 0:1:15 ForwardOrdered Regular Points,</span></span>
<span class="line"><span>  → q2di_j Sampled{Int64} 0:1:15 ForwardOrdered Regular Points</span></span>
<span class="line"><span>├───────────────────────────────────────────────────── metadata ┤</span></span>
<span class="line"><span>  Dict{String, Any} with 18 entries:</span></span>
<span class="line"><span>  &quot;dggs_radius&quot;           =&gt; 6.37101e6</span></span>
<span class="line"><span>  &quot;long_name&quot;             =&gt; &quot;monthly NDVI CMG 0.05 Deg Monthly NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_azimuth&quot; =&gt; 0</span></span>
<span class="line"><span>  &quot;scale_factor&quot;          =&gt; 0.0001</span></span>
<span class="line"><span>  &quot;dggs_rotation_lon&quot;     =&gt; 11.25</span></span>
<span class="line"><span>  &quot;dggs_id&quot;               =&gt; &quot;DGGRID ISEA4H Q2DI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polygon&quot;          =&gt; &quot;hexagon&quot;</span></span>
<span class="line"><span>  &quot;dggs_level&quot;            =&gt; 5</span></span>
<span class="line"><span>  &quot;_FillValue&quot;            =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;units&quot;                 =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;name&quot;                  =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polyhedron&quot;       =&gt; &quot;icosahedron&quot;</span></span>
<span class="line"><span>  &quot;missing_value&quot;         =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;add_offset&quot;            =&gt; 0.0</span></span>
<span class="line"><span>  &quot;dggs_projection&quot;       =&gt; &quot;isea&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_lat&quot;     =&gt; 58.2825</span></span>
<span class="line"><span>  &quot;dggs_aperture&quot;         =&gt; 4</span></span>
<span class="line"><span>  &quot;dggs_index&quot;            =&gt; &quot;Q2DI&quot;</span></span>
<span class="line"><span>├──────────────────────────────────────────────────── file size ┤ </span></span>
<span class="line"><span>  file size: 1.0 KB</span></span>
<span class="line"><span>└───────────────────────────────────────────────────────────────┘</span></span></code></pre></div><h2 id="Select-cells-and-its-neighbors" tabindex="-1">Select cells and its neighbors <a class="header-anchor" href="#Select-cells-and-its-neighbors" aria-label="Permalink to &quot;Select cells and its neighbors {#Select-cells-and-its-neighbors}&quot;">​</a></h2><p>Select a single cell using geographical coordinates (lon, lat):</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> p[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">6</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ndvi</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">a[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">11.586</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50.927</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭─────────────────────────────────────────────────╮</span></span>
<span class="line"><span>│ 216-element YAXArray{Union{Missing, Float32},1} │</span></span>
<span class="line"><span>├─────────────────────────────────────────────────┴────────────────────── dims ┐</span></span>
<span class="line"><span>  ↓ Ti Sampled{Dates.DateTime} [2001-01-01T00:00:00, …, 2018-12-01T00:00:00] ForwardOrdered Irregular Points</span></span>
<span class="line"><span>├──────────────────────────────────────────────────────────────────── metadata ┤</span></span>
<span class="line"><span>  Dict{String, Any} with 18 entries:</span></span>
<span class="line"><span>  &quot;dggs_radius&quot;           =&gt; 6.37101e6</span></span>
<span class="line"><span>  &quot;long_name&quot;             =&gt; &quot;monthly NDVI CMG 0.05 Deg Monthly NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_azimuth&quot; =&gt; 0</span></span>
<span class="line"><span>  &quot;scale_factor&quot;          =&gt; 0.0001</span></span>
<span class="line"><span>  &quot;dggs_rotation_lon&quot;     =&gt; 11.25</span></span>
<span class="line"><span>  &quot;dggs_id&quot;               =&gt; &quot;DGGRID ISEA4H Q2DI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polygon&quot;          =&gt; &quot;hexagon&quot;</span></span>
<span class="line"><span>  &quot;dggs_level&quot;            =&gt; 6</span></span>
<span class="line"><span>  &quot;_FillValue&quot;            =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;units&quot;                 =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;name&quot;                  =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polyhedron&quot;       =&gt; &quot;icosahedron&quot;</span></span>
<span class="line"><span>  &quot;missing_value&quot;         =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;add_offset&quot;            =&gt; 0.0</span></span>
<span class="line"><span>  &quot;dggs_projection&quot;       =&gt; &quot;isea&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_lat&quot;     =&gt; 58.2825</span></span>
<span class="line"><span>  &quot;dggs_aperture&quot;         =&gt; 4</span></span>
<span class="line"><span>  &quot;dggs_index&quot;            =&gt; &quot;Q2DI&quot;</span></span>
<span class="line"><span>├─────────────────────────────────────────────────────────────────── file size ┤ </span></span>
<span class="line"><span>  file size: 864.0 bytes</span></span>
<span class="line"><span>└──────────────────────────────────────────────────────────────────────────────┘</span></span></code></pre></div><p>Select the same cell using DGGS coordinates (n,i,j):</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">a[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">30</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭─────────────────────────────────────────────────╮</span></span>
<span class="line"><span>│ 216-element YAXArray{Union{Missing, Float32},1} │</span></span>
<span class="line"><span>├─────────────────────────────────────────────────┴────────────────────── dims ┐</span></span>
<span class="line"><span>  ↓ Ti Sampled{Dates.DateTime} [2001-01-01T00:00:00, …, 2018-12-01T00:00:00] ForwardOrdered Irregular Points</span></span>
<span class="line"><span>├──────────────────────────────────────────────────────────────────── metadata ┤</span></span>
<span class="line"><span>  Dict{String, Any} with 18 entries:</span></span>
<span class="line"><span>  &quot;dggs_radius&quot;           =&gt; 6.37101e6</span></span>
<span class="line"><span>  &quot;long_name&quot;             =&gt; &quot;monthly NDVI CMG 0.05 Deg Monthly NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_azimuth&quot; =&gt; 0</span></span>
<span class="line"><span>  &quot;scale_factor&quot;          =&gt; 0.0001</span></span>
<span class="line"><span>  &quot;dggs_rotation_lon&quot;     =&gt; 11.25</span></span>
<span class="line"><span>  &quot;dggs_id&quot;               =&gt; &quot;DGGRID ISEA4H Q2DI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polygon&quot;          =&gt; &quot;hexagon&quot;</span></span>
<span class="line"><span>  &quot;dggs_level&quot;            =&gt; 6</span></span>
<span class="line"><span>  &quot;_FillValue&quot;            =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;units&quot;                 =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;name&quot;                  =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polyhedron&quot;       =&gt; &quot;icosahedron&quot;</span></span>
<span class="line"><span>  &quot;missing_value&quot;         =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;add_offset&quot;            =&gt; 0.0</span></span>
<span class="line"><span>  &quot;dggs_projection&quot;       =&gt; &quot;isea&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_lat&quot;     =&gt; 58.2825</span></span>
<span class="line"><span>  &quot;dggs_aperture&quot;         =&gt; 4</span></span>
<span class="line"><span>  &quot;dggs_index&quot;            =&gt; &quot;Q2DI&quot;</span></span>
<span class="line"><span>├─────────────────────────────────────────────────────────────────── file size ┤ </span></span>
<span class="line"><span>  file size: 864.0 bytes</span></span>
<span class="line"><span>└──────────────────────────────────────────────────────────────────────────────┘</span></span></code></pre></div><p><img src="`+t+`" alt=""> A 3-ring and a 2-disk around a center cell</p><p>Select a 2-disk containing all neighboring cells that are at least k cells apart including a given center cell:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">a[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">11.586</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50.927</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭───────────────────────────────────────────╮</span></span>
<span class="line"><span>│ 7×216 YAXArray{Union{Missing, Float32},2} │</span></span>
<span class="line"><span>├───────────────────────────────────────────┴──────────────────────────── dims ┐</span></span>
<span class="line"><span>  ↓ q2di_k Sampled{Int64} 1:7 ForwardOrdered Regular Points,</span></span>
<span class="line"><span>  → Ti     Sampled{Dates.DateTime} [2001-01-01T00:00:00, …, 2018-12-01T00:00:00] ForwardOrdered Irregular Points</span></span>
<span class="line"><span>├──────────────────────────────────────────────────────────────────── metadata ┤</span></span>
<span class="line"><span>  Dict{String, Any}()</span></span>
<span class="line"><span>├─────────────────────────────────────────────────────────────────── file size ┤ </span></span>
<span class="line"><span>  file size: 5.91 KB</span></span>
<span class="line"><span>└──────────────────────────────────────────────────────────────────────────────┘</span></span></code></pre></div><p>This will introduce a new dimension <code>q2di_k</code> iterating over all neighbors. The ordering of cells within this dimension is deterministic but not further specified.</p><p>Select a 3-ring of cells having the same distance to the center cell:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">a[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">11.586</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50.927</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭────────────────────────────────────────────╮</span></span>
<span class="line"><span>│ 12×216 YAXArray{Union{Missing, Float32},2} │</span></span>
<span class="line"><span>├────────────────────────────────────────────┴─────────────────────────── dims ┐</span></span>
<span class="line"><span>  ↓ q2di_k Sampled{Int64} 1:12 ForwardOrdered Regular Points,</span></span>
<span class="line"><span>  → Ti     Sampled{Dates.DateTime} [2001-01-01T00:00:00, …, 2018-12-01T00:00:00] ForwardOrdered Irregular Points</span></span>
<span class="line"><span>├──────────────────────────────────────────────────────────────────── metadata ┤</span></span>
<span class="line"><span>  Dict{String, Any}()</span></span>
<span class="line"><span>├─────────────────────────────────────────────────────────────────── file size ┤ </span></span>
<span class="line"><span>  file size: 10.12 KB</span></span>
<span class="line"><span>└──────────────────────────────────────────────────────────────────────────────┘</span></span></code></pre></div><h2 id="long-and-short-syntax" tabindex="-1">long and short syntax <a class="header-anchor" href="#long-and-short-syntax" aria-label="Permalink to &quot;long and short syntax {#long-and-short-syntax}&quot;">​</a></h2><p>Selection on both spatial and non-spatial dimensions can be performed using keyword-based arguments on pyramids, ndvis, and arrays:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p[id</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:ndvi</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, Time</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, level</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, lon</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">11.586</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, lat</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50.927</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭───────────────────────────────────────────────────╮</span></span>
<span class="line"><span>│ 0-dimensional YAXArray{Union{Missing, Float32},0} │</span></span>
<span class="line"><span>├───────────────────────────────────────────────────┴ metadata ┐</span></span>
<span class="line"><span>  Dict{String, Any} with 18 entries:</span></span>
<span class="line"><span>  &quot;dggs_radius&quot;           =&gt; 6.37101e6</span></span>
<span class="line"><span>  &quot;long_name&quot;             =&gt; &quot;monthly NDVI CMG 0.05 Deg Monthly NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_azimuth&quot; =&gt; 0</span></span>
<span class="line"><span>  &quot;scale_factor&quot;          =&gt; 0.0001</span></span>
<span class="line"><span>  &quot;dggs_rotation_lon&quot;     =&gt; 11.25</span></span>
<span class="line"><span>  &quot;dggs_id&quot;               =&gt; &quot;DGGRID ISEA4H Q2DI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polygon&quot;          =&gt; &quot;hexagon&quot;</span></span>
<span class="line"><span>  &quot;dggs_level&quot;            =&gt; 5</span></span>
<span class="line"><span>  &quot;_FillValue&quot;            =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;units&quot;                 =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;name&quot;                  =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polyhedron&quot;       =&gt; &quot;icosahedron&quot;</span></span>
<span class="line"><span>  &quot;missing_value&quot;         =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;add_offset&quot;            =&gt; 0.0</span></span>
<span class="line"><span>  &quot;dggs_projection&quot;       =&gt; &quot;isea&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_lat&quot;     =&gt; 58.2825</span></span>
<span class="line"><span>  &quot;dggs_aperture&quot;         =&gt; 4</span></span>
<span class="line"><span>  &quot;dggs_index&quot;            =&gt; &quot;Q2DI&quot;</span></span>
<span class="line"><span>├─────────────────────────────────────────────────── file size ┤ </span></span>
<span class="line"><span>  file size: 4.0 bytes</span></span>
<span class="line"><span>└──────────────────────────────────────────────────────────────┘</span></span></code></pre></div><p>which is equivalent to:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ndvi[Time</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">11.586</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50.927</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭───────────────────────────────────────────────────╮</span></span>
<span class="line"><span>│ 0-dimensional YAXArray{Union{Missing, Float32},0} │</span></span>
<span class="line"><span>├───────────────────────────────────────────────────┴ metadata ┐</span></span>
<span class="line"><span>  Dict{String, Any} with 18 entries:</span></span>
<span class="line"><span>  &quot;dggs_radius&quot;           =&gt; 6.37101e6</span></span>
<span class="line"><span>  &quot;long_name&quot;             =&gt; &quot;monthly NDVI CMG 0.05 Deg Monthly NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_azimuth&quot; =&gt; 0</span></span>
<span class="line"><span>  &quot;scale_factor&quot;          =&gt; 0.0001</span></span>
<span class="line"><span>  &quot;dggs_rotation_lon&quot;     =&gt; 11.25</span></span>
<span class="line"><span>  &quot;dggs_id&quot;               =&gt; &quot;DGGRID ISEA4H Q2DI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polygon&quot;          =&gt; &quot;hexagon&quot;</span></span>
<span class="line"><span>  &quot;dggs_level&quot;            =&gt; 5</span></span>
<span class="line"><span>  &quot;_FillValue&quot;            =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;units&quot;                 =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;name&quot;                  =&gt; &quot;NDVI&quot;</span></span>
<span class="line"><span>  &quot;dggs_polyhedron&quot;       =&gt; &quot;icosahedron&quot;</span></span>
<span class="line"><span>  &quot;missing_value&quot;         =&gt; -9999.0</span></span>
<span class="line"><span>  &quot;add_offset&quot;            =&gt; 0.0</span></span>
<span class="line"><span>  &quot;dggs_projection&quot;       =&gt; &quot;isea&quot;</span></span>
<span class="line"><span>  &quot;dggs_rotation_lat&quot;     =&gt; 58.2825</span></span>
<span class="line"><span>  &quot;dggs_aperture&quot;         =&gt; 4</span></span>
<span class="line"><span>  &quot;dggs_index&quot;            =&gt; &quot;Q2DI&quot;</span></span>
<span class="line"><span>├─────────────────────────────────────────────────── file size ┤ </span></span>
<span class="line"><span>  file size: 4.0 bytes</span></span>
<span class="line"><span>└──────────────────────────────────────────────────────────────┘</span></span></code></pre></div>`,39)]))}const k=a(l,[["render",e]]);export{r as __pageData,k as default};
