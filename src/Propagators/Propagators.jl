"""
    Propagators

Provides orbit propagators for the two-body problem, 
and the n-body problem.
"""
module Propagators

using ..CommonTypes
using ..NBody
using ..TwoBody

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm, normalize
using ComponentArrays
using StaticArrays
using OrdinaryDiffEq

export  TwobodyPropagationResult, 
        MultibodyPropagationResult, 
        propagate

include("propagate_twobody.jl")
include("propagate_nbody.jl")

end