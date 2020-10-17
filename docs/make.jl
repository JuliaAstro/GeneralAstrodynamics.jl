push!(LOAD_PATH,"../src/")
using Documenter
using UnitfulAstrodynamics

makedocs(modules=[UnitfulAstrodynamics],
    format=Documenter.HTML(),
    sitename="UnitfulAstrodynamics.jl",
    pages=[
        "Getting Started" => "index.md",
        "Index" => Any[
            "About" => "Overview/about.md",
            "Table of Contents" => "Overview/contents.md"
        ],
        "`TwoBody`" => Any[
            "Data Structures and Types" => "TwoBody/types.md",
            "Functions" => "TwoBody/functions.md"
        ],
        "`NBody`" => Any[
            "Data Structures and Types" => "NBody/types.md",
            "Functions" => "NBody/functions.md"
        ],
        "`Propagators`" => Any[
            "Data Structures and Types" => "Propagators/types.md",
            "Functions" => "Propagators/functions.md"
        ],
        "`Plots`" => Any[
            "Functions" => "Plots/functions.md"
        ],
        "Common `AbstractTypes`" => Any[
            "Types" => "AbstractTypes/types.md"
        ]
    ]
)

deploydocs(
    target = "build",
    repo="github.com/cadojo/UnitfulAstrodynamics.jl.git",
    branch = "gh-pages",
    deps   = nothing,
    make   = nothing,
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#", "v#.#.#],
)
