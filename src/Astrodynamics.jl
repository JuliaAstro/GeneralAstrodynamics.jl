"""
    Astrodynamics

Julia package developed alongside ENAE601 at the University of Maryland.
Includes structures and functions to handle common Astrodynamics problems.
"""
module Astrodynamics

# References:
# 
# I troubleshooted Julia's module scope/export rules for a while.
# Ultimately, I referenced [1] as an example for how to structure 
# submodules, and how to using `Reexport.jl` to export required
# package dependencies.
#
# [1] https://github.com/JuliaAstro/AstroBase.jl/blob/master/src/AstroBase.jl

using Reexport

include("AbstractOrbits/AbstractOrbits.jl")
include("TwoBody/TwoBody.jl")
include("NBody/NBody.jl")

@reexport using .AbstractOrbits
@reexport using .TwoBody
@reexport using .NBody

end # module
