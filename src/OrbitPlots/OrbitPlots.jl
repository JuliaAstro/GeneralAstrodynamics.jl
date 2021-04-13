"""
Uses `Plots.jl` to provide position and velocity 
plots for R2BP, CR3BP, and NBP orbits, as well
as lagrange plots and zero velocity curves.
"""
module OrbitPlots

# Module Exports
export lagrangeplot

# Module Dependencies
using Plots
using LinearAlgebra

using ..Orbits

# Source Files
include("CR3BP/LagrangePlot.jl")
    
end