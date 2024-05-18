using Documenter
using DocStringExtensions
using GeneralAstrodynamics 
using DifferentialEquations
using Plots

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
    devbranch = "main",
    versions = ["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
