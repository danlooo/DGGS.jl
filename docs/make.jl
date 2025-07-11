using DGGS
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(DGGS, :DocTestSetup, :(using DGGS); recursive=true)

makedocs(;
    modules=[DGGS],
    authors="Daniel Loos",
    sitename="DGGS.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/danlooo/DGGS.jl",
        devbranch="main",
        devurl="dev",
        # clean_md_output=false
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "User Guide" => [
            "Get Started" => "get_started.md",
            "Background" => "background.md",
        ],
        "Tutorials" => [
            "Sentinel-2 NDVI" => "sentinel-2-ndvi.md",
            "NetCDF Climate" => "netcdf-climate.md"
        ],
    ]
)

DocumenterVitepress.deploydocs(;
    repo="github.com/danlooo/DGGS.jl",
    target="build", # this is where Vitepress stores its output
    devbranch="main",
    branch="gh-pages",
    push_preview=true
)
