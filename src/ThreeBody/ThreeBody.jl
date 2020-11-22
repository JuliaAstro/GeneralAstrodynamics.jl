"""
Handles calculations relevant to the Circular Restricted
Three Body Problem.
"""
module ThreeBody

include("../Misc/DocStringExtensions.jl")

using Reexport

@reexport using ..CommonTypes
using ..TwoBody

using LinearAlgebra: norm, cross, ×, dot, ⋅
using StaticArrays: SVector, @SVector, SMatrix, @SMatrix

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