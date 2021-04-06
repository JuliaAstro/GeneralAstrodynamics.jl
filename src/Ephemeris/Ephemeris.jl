"""
Module for loading Cartesian ephemeris data
from ASCII files, e.g. plaintext HORIZONS
outputs.
"""
module Ephemeris

export loadascii, interpolator

using Reexport
@reexport using ..CommonTypes

using DelimitedFiles, Interpolations

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

include("LoadASCII.jl")

end