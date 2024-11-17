using Documenter
using AstrodynamicalSolvers

makedocs(
    sitename = "AstrodynamicalSolvers.jl",
    format = Documenter.HTML(),
    modules = [AstrodynamicalSolvers],
    authors = "Joey Carpinelli",
    pages = [
        "Getting Started" => "index.md",
        "Reference" =>
            ["`AstrodynamicalSolvers`" => "reference.md", "`CR3BSolvers`" => "cr3bp.md"],
    ],
)


deploydocs(
    target = "build",
    repo = "github.com/JuliaAstro/AstrodynamicalSolvers.jl.git",
    branch = "docs/astrodynamical-solvers",
    devbranch = "main",
    tag_prefix = "AstrodynamicalSolvers-",
)