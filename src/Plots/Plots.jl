"""
    Plots

Uses `Plots.jl` to plot propagation results for the n-body problem, 
and the two-body problem. Uses `Plots.jl`
"""
module Plots

using ..CommonTypes
using ..TwoBody
using ..NBody
using ..Propagators

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using ComponentArrays
using StaticArrays
using Plots, Plots.PlotMeasures

export plot

include("plot_twobody.jl")
include("plot_nbody.jl")


end