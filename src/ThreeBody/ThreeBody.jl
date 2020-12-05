"""
Handles calculations relevant to the Circular Restricted
Three Body Problem.
"""
module ThreeBody

using Reexport

@reexport using ..CommonTypes
using ..TwoBody

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

using LinearAlgebra: norm, cross, ×, dot, ⋅
using StaticArrays: SVector, @SVector, SMatrix, @SMatrix
using Roots

export ThreeBodySystem
export time_scale_factor,
       nondimensionalize_length,
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
       lagrange,
       nondimensional_radius,  
       inertial, 
       synodic,
       convert,
       promote
       
include("ThreeBodyStates.jl")
include("ThreeBodyCalculations.jl")

end