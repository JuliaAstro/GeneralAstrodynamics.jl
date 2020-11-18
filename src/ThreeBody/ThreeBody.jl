"""
    ThreeBody

Handles calculations relevant to the Circular Restricted
Three Body Problem.
"""
module ThreeBody

using ..CommonTypes
using ..TwoBody

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using ComponentArrays
using StaticArrays

@reexport using Unitful, UnitfulAstro, UnitfulAngles

export ThreeBodySystem
export rendimensionalize, 
       nondimensionalize, 
       potential_energy, 
       jacobi_constant,
       position,  
       inertial, 
       synodic,
       threebody_tic!

include("ThreeBodyStates.jl")
include("ThreeBodyCalculations.jl")

end