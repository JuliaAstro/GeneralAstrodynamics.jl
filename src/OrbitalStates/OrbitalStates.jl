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

using ArrayInterface: Cartesian
export StateVector, CartesianState, KeplerianState
export OrbitalState
export lengthunit, timeunit, angularunit

using Unitful
using StaticArrays
using LabelledArrays
using ArrayInterface
using DocStringExtensions

# using ..OrbitalFrames

include("StateVectors.jl")
include("OrbitDescriptions.jl")

end # module
