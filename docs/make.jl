using DGGS
using Documenter

DocMeta.setdocmeta!(DGGS, :DocTestSetup, :(using DGGS); recursive=true)

makedocs(;
    modules=[DGGS],
    authors="Daniel Loos",
    repo="https://github.com/danlooo/DGGS.jl/blob/{commit}{path}#{line}",
    sitename="DGGS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://danlooo.github.io/DGGS.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Background" => "background.md",
        "API" => "api.md",
        "FAQ" => "faq.md"
    ]
)

deploydocs(;
    repo="github.com/danlooo/DGGS.jl",
    devbranch="main"
)
