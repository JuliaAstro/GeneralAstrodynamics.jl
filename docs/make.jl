using Documenter
using GeneralAstrodynamics 


makedocs(
    format=Documenter.HTML(),
    sitename="GeneralAstrodynamics.jl",
    authors = "Joey Carpinelli",
    pages=[
        "Quick Start" => [
            "Getting Started" => "index.md",
            "Docstrings" => "docstrings.md"
        ]
    ]
)

deploydocs(
    target = "build",
    repo="github.com/cadojo/GeneralAstrodynamics.jl.git",
    branch = "gh-pages",
    deps   = nothing,
    make   = nothing,
    devbranch = "main",
    versions = ["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
