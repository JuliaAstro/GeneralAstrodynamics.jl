using Documenter
using Quarto

Quarto.render(joinpath(@__DIR__, "src"))

deploydocs(
    target = "build",
    repo = "github.com/JuliaAstro/GeneralAstrodynamics.jl",
    branch = "docs/astrodynamical-solvers",
    devbranch = "main",
    tag_prefix = "AstrodynamicalSolvers-",
)