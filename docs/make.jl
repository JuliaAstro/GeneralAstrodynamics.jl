using Documenter
using GeneralAstrodynamics 


makedocs(
    modules=[
        GeneralAstrodynamics, 
        GeneralAstrodynamics.AstrodynamicsCore,
        GeneralAstrodynamics.Propagators,
        GeneralAstrodynamics.OrbitPlots,
        GeneralAstrodynamics.Ephemeris
    ],
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
