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

export BodyState, System, propagate, MultibodyPropagationResult

include("states.jl")
include("propagator.jl")
include("plot.jl")


end