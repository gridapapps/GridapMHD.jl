using Documenter, GridapMHD

makedocs(;
    modules=[GridapMHD],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jesusbonilla/GridapMHD.jl/blob/{commit}{path}#L{line}",
    sitename="GridapMHD.jl",
    authors="Jesús Bonilla, Large Scale Scientific Computing",
    assets=String[],
)

deploydocs(;
    repo="github.com/jesusbonilla/GridapMHD.jl",
)
