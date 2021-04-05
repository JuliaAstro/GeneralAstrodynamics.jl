"""
Provides structures for major and minor bodies in our
Solar System, and beyond!
"""
module CelestialBodies

using Reexport
@reexport using ..CommonTypes

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

include("CelestialTypes.jl")
include("Constants.jl")

end