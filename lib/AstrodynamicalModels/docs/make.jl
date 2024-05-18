using Latexify
using Documenter
using ModelingToolkit
using AstrodynamicalModels

makedocs(
    sitename="AstrodynamicalModels",
    format=Documenter.HTML(),
    modules=[AstrodynamicalModels],
    pages=[
        "Overview" => [
            "Getting Started" => "index.md",
            "Docstrings" => "docstrings.md"
        ],
        "Models" => [
            "R2BP" => "R2BP.md",
            "CR3BP" => "CR3BP.md",
            "NBP" => "NBP.md",
            "Entry" => "Entry.md",
            "Attitude" => "Attitude.md",
        ]
    ]
)

deploydocs(
    target="build",
    repo="github.com/cadojo/AstrodynamicalModels.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions=["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
