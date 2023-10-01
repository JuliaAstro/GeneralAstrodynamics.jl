using Documenter
using AstrodynamicalSolvers

makedocs(
    sitename="AstrodynamicalSolvers.jl",
    format=Documenter.HTML(),
    modules=[AstrodynamicalSolvers],
    authors="Joey Carpinelli",
    pages=[
        "Getting Started" => "index.md",
        "Reference" => [
            "`AstrodynamicalSolvers`" => "reference.md",
            "`CR3BSolvers`" => "cr3bp.md",
        ]
    ]
)

deploydocs(
    target="build",
    repo="github.com/cadojo/AstrodynamicalSolvers.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions=["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
