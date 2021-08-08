"""
A module which provides types for common 
state representations in astrodynamics, 
including Cartesian and Keplerian states,
and orbital states which include epoch and 
coordinate frame information.

# Extended Help

**Exports**

$(EXPORTS)

**Imports**

$(IMPORTS)
"""
module States

export CartesianState, CartesianStateWithSTM, KeplerianState
export R2BPParameters, CR3BPParameters
export Orbit, state, system, epoch
export position, velocity
export R2BPOrbit, KeplerianR2BPOrbit, CartesianR2BPOrbit, CartesianOrbitWithSTM, CR3BPOrbit
export lengthunit, timeunit, angularunit, velocityunit, massparamunit, name
export Sun, Mercury, Venus, Earth, Moon, Luna, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
export SunVenus, SunEarth, EarthMoon, SunMars, SunJupiter, SunSaturn
export model, vectorfield

import Dates: now
import AstroTime: TAIEpoch
import Requires: @require
import LinearAlgebra: norm

using Unitful, UnitfulAstro
using StaticArrays
using ArrayInterface
using LabelledArrays
using DocStringExtensions

using ..CoordinateFrames
using ..Calculations

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

include(joinpath("Common","ParameterizedLabelledArrays.jl"))
include(joinpath("States", "StateVectors.jl"))
include(joinpath("Systems", "ParameterVectors.jl"))
include(joinpath("Orbits", "OrbitDescriptions.jl"))
include(joinpath("Systems", "R2BPSystems.jl"))
include(joinpath("Systems", "CR3BPSystems.jl"))

include("Calculations/Calculations.jl")

function __init__()
    @require AstrodynamicalModels="4282b555-f590-4262-b575-3e516e1493a7" include(joinpath(@__DIR__, "Hooks", "AstrodynamicalModels.jl"))
    @require SymbolicUtils="d1185830-fcd6-423d-90d6-eec64667417b" include(joinpath(@__DIR__, "Hooks", "SymbolicUtils.jl"))
end

end # module
