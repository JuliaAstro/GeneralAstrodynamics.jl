"""
A module which provides types for common 
state representations in astrodynamics, 
including Cartesian and Keplerian states,
and orbital states which include epoch and 
coordinate frame information.

# Extended help

**Exports**

$(EXPORTS)

**Imports**

$(IMPORTS)
"""
module OrbitalStates

using StaticArrays: getproperty
export StateVector, CartesianState, KeplerianState
export ParameterVector, R2BPParameters, CR3BPParameters
export AbstractOrbit, Orbit, state, system
export R2BPOrbit, CR3BPOrbit
export lengthunit, timeunit, angularunit

import Dates: now
using AstroTime
using Unitful
using Requires
using StaticArrays
using LabelledArrays
using ArrayInterface
using DocStringExtensions
import ArrayInterface: Cartesian

function __init__()
    @require AstrodynamicalModels="4282b555-f590-4262-b575-3e516e1493a7" include(joinpath(@__DIR__, "Hooks", "AstrodynamicalModels.jl"))
end

using ..OrbitalFrames

include("ParameterizedLabelledArrays.jl")
include("StateVectors.jl")
include("ParameterVectors.jl")
include("OrbitDescriptions.jl")

end # module
