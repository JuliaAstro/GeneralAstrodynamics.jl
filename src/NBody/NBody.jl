"""
Handles the non-relativistic NBody problem for planets,
and other celestial bodies.
"""
module NBody

using Reexport
@reexport using ..AstrodynamicsCore

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

using StaticArrays: SVector, @SVector, SMatrix, @SMatrix

export Body, NBodySystem, system_energy, 
       system_angular_momentum, promote, convert,
       Float16, Float32, Float64, BigFloat,
       length, getindex

include("NBodyStates.jl")
include("NBodyCalculations.jl")


end