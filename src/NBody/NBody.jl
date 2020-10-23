"""
    NBody

Handles the non-relativistic NBody problem for planets,
and other celestial bodies.
"""
module NBody

using ..CommonTypes

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using ComponentArrays
using StaticArrays

@reexport using Unitful, UnitfulAstro, UnitfulAngles

export Body, MultibodySystem, system_energy, 
       system_angular_momentum

include("multibody_states.jl")
include("multibody_calculations.jl")


end