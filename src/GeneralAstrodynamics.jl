"""
Structures and functions for handling common Astrodynamics problems! 🚀
"""
module GeneralAstrodynamics

using Reexport 

include("AstrodynamicsCore/AstrodynamicsCore.jl")
include("Propagators/Propagators.jl")
include("OrbitPlots/OrbitPlots.jl")
include("Ephemeris/Ephemeris.jl")

@reexport using .AstrodynamicsCore
@reexport using .Propagators
@reexport using .OrbitPlots
@reexport using .Ephemeris

end # module
