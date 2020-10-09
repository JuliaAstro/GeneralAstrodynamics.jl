"""
    Plots

Uses `Plots.jl` to plot propagation results for the n-body problem, 
and the two-body problem. Uses `Plots.jl`
"""
module Plots

using ..AbstractTypes
using ..TwoBody
using ..NBody
using ..Propagators

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using DifferentialEquations
using ComponentArrays
using StaticArrays

@reexport using Unitful, UnitfulAstro, UnitfulAngles
@reexport using Plots, Plots.PlotMeasures

export plot3d

include("plot_twobody.jl")
include("plot_nbody.jl")


end