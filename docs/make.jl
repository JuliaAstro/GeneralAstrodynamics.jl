using Documenter
using Orbits 


makedocs(
    modules=[
       Orbits 
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
    repo="github.com/cadojo/Orbits.jl.git",
    branch = "gh-pages",
    deps   = nothing,
    make   = nothing,
    devbranch = "main",
    versions = ["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
