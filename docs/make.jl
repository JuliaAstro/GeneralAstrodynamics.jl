using Documenter
using Quarto
using Git

Quarto.render(joinpath(@__DIR__, "src"))

function clone(pkg; branch = "docs/$pkg")
    run(Git.git(["clone", "--branch", branch, "--depth", "1", "https://github.com/JuliaAstro/GeneralAstrodynamics.jl", joinpath(@__DIR__, "build", "lib", "$pkg.jl")]))
    rm(joinpath("build", "lib", "$pkg.jl", ".git"), force=true, recursive=true)
end

packages = [
    name
    for name in readdir(joinpath(@__DIR__, "..", "lib"))
    if isdir(joinpath(@__DIR__, "..", "lib", name))
]

@info "packages = $packages"
for package in packages
    rm(joinpath(@__DIR__, "build", "lib", "$package.jl"); force=true, recursive=true)
    clone(package)
end

Documenter.deploydocs(repo = "github.com/JuliaAstro/GeneralAstrodynamics.jl", push_preview=true)
