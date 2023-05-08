using DGGS
using Documenter

DocMeta.setdocmeta!(DGGS, :DocTestSetup, :(using DGGS); recursive=true)

makedocs(;
    modules=[DGGS],
    authors="Daniel Loos",
    repo="https://github.com/dloos/DGGS.jl/blob/{commit}{path}#{line}",
    sitename="DGGS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dloos.github.io/DGGS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dloos/DGGS.jl",
    devbranch="main",
)
