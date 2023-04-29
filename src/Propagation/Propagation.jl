"""
A module which provides wrappers around
`DifferentialEquations` solvers and
`AstrodynamicalModels` for orbit
propagation, iterative Halo orbit solvers,
and manifold calculations.

# Extended Help

**Exports**

$(EXPORTS)

**Imports**

$(IMPORTS)
"""
module Propagation

export Trajectory, ODEProblem, propagate
export Manifold, EnsembleProblem
export initialstate, initialepoch, solution
export isperiodic, halo, monodromy
export manifold, perturb
export stable_eigenvector, unstable_eigenvector

using ..States
using ..Calculations
using ..CoordinateFrames

import Dates: now
import LinearAlgebra: I
import Roots: find_zero

using Unitful
using Requires
using AstroTime
using Distributed
using StaticArrays
using LinearAlgebra
using ModelingToolkit
using DocStringExtensions
using AstrodynamicalModels
using DifferentialEquations

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)
                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)
                               $(DOCSTRING)
                               """

include("Trajectories.jl")
include("Propagators.jl")
include("Halos.jl")
include("Manifolds.jl")

function __init__()
    @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" begin
        DataFrames.DataFrame(traj::Trajectory) = DataFrames.DataFrame(solution(traj))
    end
end

end # module
