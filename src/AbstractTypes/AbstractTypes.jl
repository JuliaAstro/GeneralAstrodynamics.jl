"""
    AbstractTypes

Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module AbstractTypes

using Reexport
@reexport using Unitful, UnitfulAstro

export OrbitalSystem, OrbitalState, PropagationResult

"""
    AbstractSystem

Abstract type describing all states in select Astrodynamics problems.
"""
abstract type OrbitalSystem end

"""
    AbstractState

Abstract type for orbital states.
"""
abstract type OrbitalState end

"""
    PropagationResult

Abstract type describing a collection of states resulting from 
"""
abstract type PropagationResult end

end