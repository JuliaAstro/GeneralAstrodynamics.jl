"""
    Propagators

Provides orbit propagators for the two-body problem, 
and the n-body problem.
"""
module Propagators

using ..AbstractTypes
using ..NBody
using ..TwoBody

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using DifferentialEquations
using ComponentArrays
using StaticArrays

@reexport using Unitful, UnitfulAstro, UnitfulAngles

export  TwobodyPropagationResult, 
        MultibodyPropagationResult, 
        propagate

include("propagate_twobody.jl")
include("propagate_nbody.jl")

end