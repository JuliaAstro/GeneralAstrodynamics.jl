"""
Provides orbit propagators for the two-body problem, 
and the n-body problem.
"""
module Propagators

using Reexport 
@reexport using ..CommonTypes

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

using ..NBody
using ..TwoBody
using ..ThreeBody 

using OrdinaryDiffEq
using LinearAlgebra: norm, normalize, cross, ×, dot, ⋅
using ComponentArrays

export  TwobodyPropagationResult, 
        ThreeBodyPropagationResult,
        MultibodyPropagationResult, 
        propagate,
        twobody_tic!,
        threebody_tic!,
        nbody_tic,
        show

include("PropagateTwoBody.jl")
include("PropagateThreeBody.jl")
include("PropagateNBody.jl")

end