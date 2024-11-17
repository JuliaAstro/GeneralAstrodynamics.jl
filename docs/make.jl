using Documenter
using DocStringExtensions
using GeneralAstrodynamics
using DifferentialEquations
using Plots

makedocs(
    format = Documenter.HTML(size_threshold = nothing),
    sitename = "GeneralAstrodynamics.jl",
    authors = "Joey Carpinelli",
    pages = [
        "Quick Start" => ["Getting Started" => "index.md", "Docstrings" => "docstrings.md"],
    ],
)

deploydocs(
    target = "build",
    repo = "github.com/JuliaAstro/GeneralAstrodynamics.jl.git",
    branch = "docs/general-astrodynamics",
    devbranch = "main",
)


using MultiDocumenter

clonedir = mktempdir()

content = [
    MultiDocumenter.MultiDocRef(
        upstream = joinpath(clonedir, "GeneralAstrodynamics.jl"),
        path = "docs",
        name = "GeneralAstrodynamics.jl",
        branch = "docs/general-astrodynamics",
        giturl = "https://github.com/JuliaAstro/GeneralAstrodynamics.jl.git",
        fix_canonical_url = false,
    ),
    MultiDocumenter.MultiDocRef(
        upstream = joinpath(clonedir, "AstrodynamicalCalculations.jl"),
        path = joinpath("docs", "lib", "AstrodynamicalCalculations"),
        name = "Calculations",
        branch = "docs/astrodynamical-calculations",
        giturl = "https://github.com/JuliaAstro/GeneralAstrodynamics.jl.git",
        fix_canonical_url = false,
    ),
    MultiDocumenter.MultiDocRef(
        upstream = joinpath(clonedir, "AstrodynamicalModels.jl"),
        path = joinpath("docs", "lib", "AstrodynamicalModels"),
        name = "Models",
        branch = "docs/astrodynamical-models",
        giturl = "https://github.com/JuliaAstro/GeneralAstrodynamics.jl.git",
        fix_canonical_url = false,
    ),
    MultiDocumenter.MultiDocRef(
        upstream = joinpath(clonedir, "AstrodynamicalSolvers.jl"),
        path = joinpath("docs", "lib", "AstrodynamicalSolvers"),
        name = "Solvers",
        branch = "docs/astrodynamical-solvers",
        giturl = "https://github.com/JuliaAstro/GeneralAstrodynamics.jl.git",
        fix_canonical_url = false,
    ),
]

outpath = joinpath(@__DIR__, "build")

MultiDocumenter.make(
    outpath,
    content;
    prettyurls = true,
    search_engine = MultiDocumenter.SearchConfig(
        index_versions = ["stable", "dev"],
        engine = MultiDocumenter.FlexSearch,
    ),
    brand_image = MultiDocumenter.BrandImage(
        "https://juliaastro.org",
        "http://juliaastro.org/dev/assets/logo.svg",
    ),
)

Documenter.deploydocs(
    target = outpath,
    versions = nothing,
    repo = "github.com/JuliaAstro/GeneralAstrodynamics.jl",
)
