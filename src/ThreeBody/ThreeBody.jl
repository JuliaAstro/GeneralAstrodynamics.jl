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
export nondimensionalize_length,
       nondimensionalize_velocity,
       nondimensionalize_time,
       nondimensionalize_mass_parameter,
       nondimensionalize,
       redimensionalize_length,
       redimensionalize_velocity,
       redimensionalize_time,
       redimensionalize,
       potential_energy, 
       jacobi_constant,
       nondimensional_radius,  
       inertial, 
       synodic,
       convert,
       promote
       
include("ThreeBodyStates.jl")
include("ThreeBodyCalculations.jl")

end