using Documenter
using AstrodynamicalCalculations

makedocs(
    sitename="AstrodynamicalCalculations.jl",
    format=Documenter.HTML(),
    modules=[AstrodynamicalCalculations],
    authors = "Joey Carpinelli",
    pages=[
        "Getting Started" => "index.md",
        "R2BP Equations" => "r2bp.md",
        "CR3BP Equations" => "cr3bp.md",
    ]
)

deploydocs(
    target="build",
    repo="github.com/cadojo/AstrodynamicalCalculations.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions=["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
