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

using DifferentialEquations 
using LinearAlgebra: norm, normalize, cross, ×, dot, ⋅
using ComponentArrays
using StaticArrays: SVector, @SVector, SMatrix, @SMatrix

export  Trajectory,
        propagate,
        RestrictedTwoBodyTic!,
        RestrictedBiasedTwoBodyTic!,
        RestrictedThreeBodyTic!,
        NBodyTic!,
        show

include("Trajectory.jl")
include("PropagateTwoBody.jl")
include("PropagateThreeBody.jl")
include("PropagateNBody.jl")

end
