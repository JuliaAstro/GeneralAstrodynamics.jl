"""
Uses `Plots.jl` to plot propagation results for the n-body problem, 
and the two-body problem. Uses `Plots.jl`
"""
module AstroPlots

using Reexport
@reexport using ..CommonTypes

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

using ..TwoBody
using ..ThreeBody
using ..NBody
using ..Propagators

using Plots
using Plots.PlotMeasures

export plot, lagrangeplot

include("plot_twobody.jl")
include("plot_nbody.jl")
include("PlotThreeBody.jl")

end
