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

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles
@reexport using Plots, Plots.PlotMeasures

export Body, MultibodySystem, propagate_multibody, multibody_plot3d, 
       MultibodyPropagationResult, multibody_plot3d, system_energy, 
       system_angular_momentum

include("multibody_states.jl")
include("multibody_calculations.jl")
include("propagate_multibody.jl")
include("plot_multibody.jl")


end