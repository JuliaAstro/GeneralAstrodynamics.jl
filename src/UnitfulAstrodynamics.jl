"""
Julia package developed alongside ENAE601 at the University of Maryland.
Includes structures and functions to handle common Astrodynamics problems.
"""
module UnitfulAstrodynamics

# References:
# 
# I troubleshooted Julia's module scope/export rules for a while.
# Ultimately, I referenced [1] as an example for how to structure 
# submodules, and how to using `Reexport.jl` to export required
# package dependencies.
#
# [1] https://github.com/JuliaAstro/AstroBase.jl/blob/master/src/AstroBase.jl

using Reexport

include("Foundation/Foundation.jl")

@reexport using .Foundation

#=

include("AstrodynamicsCore/AstrodynamicsCore.jl")
include("TwoBody/TwoBody.jl")
include("ThreeBody/ThreeBody.jl")
include("NBody/NBody.jl")
include("Propagators/Propagators.jl")
include("Maneuvers/Maneuvers.jl")
include("OrbitPlots/OrbitPlots.jl")
include("Ephemeris/Ephemeris.jl")

@reexport using .AstrodynamicsCore
@reexport using .TwoBody
@reexport using .ThreeBody
@reexport using .NBody
@reexport using .Propagators
@reexport using .Maneuvers
@reexport using .OrbitPlots
@reexport using .Ephemeris

include("Misc/DocStringExtensions.jl")
include("Misc/UnitfulAliases.jl")

=#
end # module
