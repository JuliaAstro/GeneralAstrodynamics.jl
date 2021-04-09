using Documenter
using UnitfulAstrodynamics


makedocs(
    modules=[
        UnitfulAstrodynamics
    ],
    format=Documenter.HTML(),
    sitename="UnitfulAstrodynamics.jl",
    authors = "Joey Carpinelli",
    pages=[
        "Guide" => "index.md",
        "Overview" => Any[
            "About" => "Overview/about.md",
            "Getting Stated" => "Overview/getting-started.md"
        ],
        "Two-body" => Any[
            "Data Structures and Types" => "TwoBody/types.md",
            "Functions" => "TwoBody/functions.md"
        ],
        "Three-body" => Any[
            "Data Structures and Types" => "ThreeBody/types.md",
            "Functions" => "ThreeBody/functions.md"
        ],
        "N-body" => Any[
            "Data Structures and Types" => "NBody/types.md",
            "Functions" => "NBody/functions.md"
        ],
        "Propagators" => Any[
            "Data Structures and Types" => "Propagators/types.md",
            "Functions" => "Propagators/functions.md"
        ],
        "Maneuvers" => Any[
            "Data Structures and Types" => "Maneuvers/types.md",
            "Functions" => "Maneuvers/functions.md"
        ],
        "Plotting" => Any[
            "Functions" => "OrbitPlots/functions.md"
        ],
        "Common Types" => Any[
            "Types" => "CommonTypes/types.md"
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
    versions = ["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
