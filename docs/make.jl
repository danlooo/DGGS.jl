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
        "Guide" => [
            "Get Started" => "get_started.md",
            "Background" => "background.md",
            "Convert" => "convert.md",
            "Select" => "select.md",
            "Plot" => "plot.md",
        ],
        "API" => "api.md",
    ]
)

DocumenterVitepress.deploydocs(;
    repo="github.com/danlooo/DGGS.jl",
    target="build", # this is where Vitepress stores its output
    devbranch="main",
    branch = "gh-pages",
    push_preview = true
)
