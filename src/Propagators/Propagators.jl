"""
Contains all functions and structures that require `DifferentialEquations` solvers.
This module provides orbit propagation, and Halo orbit solvers.
"""
module Propagators

# Module Exports
export propagate, halo, monodromy, isperiodic, potential_energy_hessian
export stable_eigenvector, unstable_eigenvector, manifold
export stable_manifold, unstable_manifold
export ODEProblem
export R2BPTic!, CR3BPTic!, CR3BPSTMTic!

# Module Dependencies
using Reexport

using Distributed
using StaticArrays
using LinearAlgebra
using SymbolicUtils
using ComponentArrays
using DifferentialEquations
using Unitful, UnitfulAstro, UnitfulAngles

using ..AstrodynamicsCore

include("R2BP/R2BPPropagators.jl")
include("CR3BP/CR3BPPropagators.jl")
include("CR3BP/Halos.jl")

end