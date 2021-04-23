using Documenter
using SimpleAstrodynamics


makedocs(
    modules=[
        SimpleAstrodynamics
    ],
    format=Documenter.HTML(),
    sitename="SimpleAstrodynamics.jl",
    authors = "Joey Carpinelli",
    pages=[
        "Home" => "index.md"
    ]
)

deploydocs(
    target = "build",
    repo="github.com/cadojo/SimpleAstrodynamics.jl.git",
    branch = "gh-pages",
    deps   = nothing,
    make   = nothing,
    devbranch = "main",
    versions = ["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
