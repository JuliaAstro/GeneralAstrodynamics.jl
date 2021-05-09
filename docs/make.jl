using Documenter
using GeneralAstrodynamics 


makedocs(
    modules=[
       GeneralAstrodynamics 
    ],
    format=Documenter.HTML(),
    sitename="GeneralAstrodynamics.jl",
    authors = "Joey Carpinelli",
    pages=[
        "Home" => "index.md"
    ]
)

deploydocs(
    target = "build",
    repo="github.com/cadojo/GeneralAstrodynamics.jl.git",
    branch = "gh-pages",
    deps   = nothing,
    make   = nothing,
    devbranch = "main",
    versions = ["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
