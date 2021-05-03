"""
Uses `Plots.jl` to provide position and velocity 
plots for R2BP, CR3BP, and NBP orbits, as well
as lagrange plots and zero velocity curves.
"""
module OrbitPlots

# Module Exports
export plotpositions, plotpositions!
export plotvelocities, plotvelocities!
export lagrangeplot, lagrangeplot!
export plotbodies, plotbodies!
export zerovelocityplot, zerovelocityplot!

# Module Dependencies
using Plots
using LinearAlgebra

using ..OrbitsBase

# Source Files
include("R2BP/R2BPPlots.jl")
include("CR3BP/CR3BPPlots.jl")
include("CR3BP/LagrangePlot.jl")
include("CR3BP/ZeroVelocityPlot.jl")
end