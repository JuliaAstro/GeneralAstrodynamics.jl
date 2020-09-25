"""
    NBody

Handles the non-relativistic NBody problem for planets,
and other celestial bodies.
"""
module NBody

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using DifferentialEquations
using ComponentArrays
using Plots

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles

export BodyState, System, npropagate, NBodyPropagationResult, PropagationResult

include("States.jl")
include("Propagator.jl")
include("plot.jl")


end