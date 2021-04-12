"""
Contains all functions and structures that require `DifferentialEquations` solvers.
This module provides orbit propagation, and Halo orbit solvers.
"""
module Propagators

# Module Exports
export propagate, halo, monodromy, isperiodic, potential_energy_hessian
export ODEProblem
export R2BPTic!, CR3BPTic!, CR3BPSTMTic!

# Module Dependencies
using Reexport

using StaticArrays
using LinearAlgebra
using ComponentArrays
using DifferentialEquations

using ..Orbits

include("R2BPPropagators.jl")
include("CR3BPPropagators.jl")

end