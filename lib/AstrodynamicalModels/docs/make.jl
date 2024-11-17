using Latexify
using Documenter
using ModelingToolkit
using AstrodynamicalModels

makedocs(
    sitename = "AstrodynamicalModels",
    format = Documenter.HTML(),
    modules = [AstrodynamicalModels],
    pages = [
        "Overview" => ["Getting Started" => "index.md", "Docstrings" => "docstrings.md"],
        "Models" => [
            "R2BP" => "R2BP.md",
            "CR3BP" => "CR3BP.md",
            "NBP" => "NBP.md",
            "Entry" => "Entry.md",
            "Attitude" => "Attitude.md",
        ],
    ],
)

deploydocs(
    target = "build",
    repo = "github.com/JuliaAstro/GeneralAstrodynamics.jl.git",
    branch = "docs/astrodynamical-models",
    devbranch = "main",
    tag_prefix = "AstrodynamicalModels-",
)