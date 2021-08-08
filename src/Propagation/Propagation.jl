module Propagation

export Trajectory, ODEProblem, propagate
export initialstate, initialepoch, solution

using ..States
using ..Calculations
using ..CoordinateFrames

using Unitful
using AstroTime
using DocStringExtensions
using AstrodynamicalModels
using DifferentialEquations

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    $(METHODLIST)
    """

@template (TYPES, CONSTANTS) =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

include("Trajectories.jl")
include("Propagators.jl")

end # module
