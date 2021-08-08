module Propagation

export Trajectory 

using ..States
using ..Calculations
using ..CoordinateFrames

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

include(joinpath("Trajectories", "Trajectories.jl"))
include(joinpath("R2BP", "R2BPPropagators.jl"))
include(joinpath("CR3BP", "CR3BPPropagators.jl"))

end # module
