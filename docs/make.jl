using Documenter
using AstrodynamicalCalculations

makedocs(
    sitename="AstrodynamicalCalculations",
    format=Documenter.HTML(),
    modules=[AstrodynamicalCalculations]
)

deploydocs(
    target="build",
    repo="github.com/cadojo/AstrodynamicalCalculations.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions=["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
