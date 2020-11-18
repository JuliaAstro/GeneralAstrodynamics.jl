"""
    Propagators

Provides orbit propagators for the two-body problem, 
and the n-body problem.
"""
module Propagators

using ..CommonTypes
using ..NBody
using ..TwoBody
using ..ThreeBody 

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm, normalize
using ComponentArrays
using StaticArrays
using OrdinaryDiffEq

export  TwobodyPropagationResult, 
        ThreeBodyPropagationResult,
        MultibodyPropagationResult, 
        propagate,
        twobody_tic!,
        threebody_tic!,
        nbody_tic

include("PropagateTwoBody.jl")
include("PropagateThreeBody.jl")
include("PropagateNBody.jl")

end