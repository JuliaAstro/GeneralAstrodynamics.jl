"""
    NBody

Handles the non-relativistic NBody problem for planets,
and other celestial bodies.
"""
module NBody

using ..AbstractTypes

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using DifferentialEquations
using ComponentArrays
using Plots

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles

export MultibodyState, MultibodySystem, propagate_multibody, MultibodyPropagationResult

include("multibody_states.jl")
include("propagate_multibody.jl")
include("plot_multibody.jl")


end