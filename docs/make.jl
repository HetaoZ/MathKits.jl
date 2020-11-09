using MathKits
using Documenter

makedocs(;
    modules=[MathKits],
    authors="Hetao Z.",
    repo="https://github.com/HetaoZ/MathKits.jl/blob/{commit}{path}#L{line}",
    sitename="MathKits.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HetaoZ.github.io/MathKits.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HetaoZ/MathKits.jl",
)
