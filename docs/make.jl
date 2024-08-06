using DGGS
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(DGGS, :DocTestSetup, :(using DGGS); recursive=true)

makedocs(;
    modules=[DGGS],
    authors="Daniel Loos",
    repo="https://github.com/danlooo/DGGS.jl/blob/{commit}{path}#{line}",
    sitename="DGGS.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/danlooo/DGGS.jl",
        devbranch="main",
        devurl="dev",
        clean_md_output=false
    ),
    pages=[
        "Home" => "index.md",
        "Get started" => "get_started.md",
        "Background" => "background.md",
        "API" => "api.md",
        "FAQ" => "faq.md"
    ]
)

deploydocs(;
    repo="github.com/danlooo/DGGS.jl",
    devbranch="main"
)
