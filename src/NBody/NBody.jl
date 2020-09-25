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
using RecursiveArrayTools

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles

export BodyState, System, npropagate

include("States.jl")
include("Propagator.jl")

end