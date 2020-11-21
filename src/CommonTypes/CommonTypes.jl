"""
Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module CommonTypes

using Reexport

@reexport using Unitful, UnitfulAstro

export AbstractBody, OrbitalSystem, PropagationResult

""" 
Abstract type for bodies in space: both `CelestialBody`s (in
`TwoBody.jl`), and `Body`s (in `NBody.jl`).
"""
abstract type AbstractBody end

"""
Abstract type describing all states in select Astrodynamics problems.
"""
abstract type OrbitalSystem end

"""
Abstract type describing a collection of states resulting from 
"""
abstract type PropagationResult end

end
