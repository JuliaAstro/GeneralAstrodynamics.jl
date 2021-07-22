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
module AstrodynamicalStates

export StateVector, CartesianState, KeplerianState
export ParameterVector, R2BPParameters, CR3BPParameters
export Orbit, state, system, epoch
export R2BPOrbit, KeplerianR2BPOrbit, CartesianR2BPOrbit, CR3BPOrbit
export lengthunit, timeunit, angularunit, massparamunit
export position, velocity, distance, speed
export eccentricity, semimajor_axis, inclination
export RAAN, argument_of_periapsis, true_anomoly
export massparameter, massparameters, normalized_massparameter
export primary_massparameter, secondary_massparameter

import Dates: now
import AstroTime: UTCEpoch
import Requires: @require
import LinearAlgebra: norm

using Unitful
using StaticArrays
using ArrayInterface
using LabelledArrays
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    """

@template (TYPES, CONSTANTS) =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

using ..AstrodynamicalFrames

include(joinpath("Common","ParameterizedLabelledArrays.jl"))
include(joinpath("States", "StateVectors.jl"))
include(joinpath("Systems", "ParameterVectors.jl"))
include(joinpath("Orbits", "OrbitDescriptions.jl"))

function __init__()
    @require AstrodynamicalModels="4282b555-f590-4262-b575-3e516e1493a7" include(joinpath(@__DIR__, "Hooks", "AstrodynamicalModels.jl"))
    @require SymbolicUtils="d1185830-fcd6-423d-90d6-eec64667417b" include(joinpath(@__DIR__, "Hooks", "SymbolicUtils.jl"))
end

end # module
