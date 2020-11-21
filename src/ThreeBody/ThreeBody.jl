"""
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
export redimensionalize, 
       nondimensionalize, 
       potential_energy, 
       jacobi_constant,
       position,  
       inertial, 
       synodic
       
include("ThreeBodyStates.jl")
include("ThreeBodyCalculations.jl")

end