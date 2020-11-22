"""
Handles the non-relativistic NBody problem for planets,
and other celestial bodies.
"""
module NBody

include("../Misc/DocStringExtensions.jl")

using Reexport
@reexport using ..CommonTypes

using StaticArrays: SVector, @SVector, SMatrix, @SMatrix

export Body, MultibodySystem, system_energy, 
       system_angular_momentum, promote, convert

include("multibody_states.jl")
include("multibody_calculations.jl")


end