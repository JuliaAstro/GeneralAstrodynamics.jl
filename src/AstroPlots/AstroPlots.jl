"""
Uses `Plots.jl` to plot propagation results for the n-body problem, 
and the two-body problem. Uses `Plots.jl`
"""
module AstroPlots

include("../Misc/DocStringExtensions.jl")

using Reexport
@reexport using ..CommonTypes

using ..TwoBody
using ..NBody
using ..Propagators

using Plots
using Plots.PlotMeasures

export plot

include("plot_twobody.jl")
include("plot_nbody.jl")

end
