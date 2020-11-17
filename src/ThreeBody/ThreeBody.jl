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

include("threebody_states.jl")
include("threebody_calculations.jl")


end