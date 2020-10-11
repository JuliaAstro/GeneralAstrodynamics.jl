"""
    UnitfulAstrodynamics

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

include("AbstractTypes/AbstractTypes.jl")
include("TwoBody/TwoBody.jl")
include("NBody/NBody.jl")
include("Propagators/Propagators.jl")
include("Plots/Plots.jl")

@reexport using .AbstractTypes
@reexport using .TwoBody
@reexport using .NBody
@reexport using .Propagators
@reexport using .Plots

end # module
