"""
    ThreeBody

Handles calculations relevant to the Circular Restricted
Three Body Problem.
"""
module ThreeBody

using ..CommonTypes

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using ComponentArrays
using StaticArrays

@reexport using Unitful, UnitfulAstro, UnitfulAngles

export ThreeBodySystem
export nondimensionalize, position, potential_energy, jacobi_constant

include("ThreeBodyStates.jl")
include("ThreeBodyCalculations.jl")

end