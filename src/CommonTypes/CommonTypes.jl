"""
Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module CommonTypes

include("../Misc/DocStringExtensions.jl")

using Reexport
@reexport using Unitful, UnitfulAngles, UnitfulAstro

export AbstractBody, OrbitalSystem, PropagationResult, 
       Length, Velocity, Time, Mass, MassParameter

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

@derived_dimension MassParameter Unitful.𝐋^3/Unitful.𝐓^2
"""
Custom `Unitful` dimension for gravitational parameters.
"""
const MassParameter = MassParameter

"""
Aliases for the `Unitful` length dimension.
"""
const Length = Unitful.Length

"""
Aliases for the `Unitful` velocity dimension.
"""
const Velocity = Unitful.Velocity

"""
Aliases for the `Unitful` time dimension.
"""
const Time = Unitful.Time

"""
Aliases for the `Unitful` mass dimension.
"""
const Mass = Unitful.Mass 

end
