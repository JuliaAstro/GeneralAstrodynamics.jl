using Documenter
using Quarto
using Git

Quarto.render(joinpath(@__DIR__, "src"))

function clone(pkg; branch = "docs/$pkg")
    run(Git.git(["clone", "--branch", branch, "--depth", 1, "https://github.com/JuliaAstro/GeneralAstrodynamics.jl", joinpath("build", "lib", "$pkg.jl")]))
    rm(joinpath("build", "lib", "$pkg.jl", ".git"), force=true, recursive=true)
end

packages = (
    name
    for name in readdir(joinpath(@__DIR__, "..", "lib"))
    if isdir(name)
)

for package in packages
    clone(package)
end

Documenter.deploydocs(repo = "github.com/JuliaAstro/GeneralAstrodynamics.jl", push_preview=true)
