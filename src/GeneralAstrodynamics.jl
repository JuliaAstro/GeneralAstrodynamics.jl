"""
Structures and functions for handling common Astrodynamics problems! 🚀
"""
module GeneralAstrodynamics

using Reexport 

include("AstrodynamicsCore/AstrodynamicsCore.jl")
include("Propagators/Propagators.jl")
include("OrbitPlots/OrbitPlots.jl")
include("Ephemeris/Ephemeris.jl")

include("AstrodynamicalFrames/AstrodynamicalFrames.jl")
include("AstrodynamicalStates/AstrodynamicalStates.jl")

using .AstrodynamicsCore
@reexport using .Propagators
@reexport using .OrbitPlots
@reexport using .Ephemeris

@reexport using .OrbitalFrames
@reexport using .OrbitalStates

end # module
