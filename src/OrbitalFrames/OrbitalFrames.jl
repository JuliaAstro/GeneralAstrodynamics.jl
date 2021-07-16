"""
A module which provides types, and transformations 
for provided and user-provided orbital frames, e.g.
Earth-Centered-Inertial, Heliocentric-Inertial,
Synodic.

# Extended help

**Exports**

$(EXPORTS)

**Imports**

$(IMPORTS)
"""
module OrbitalFrames

# Frames
export AstrodynamicsFrame, Inertial, Rotating
export BodycentricInertial, BarycentricInertial
export BodycentricRotating, BarycentricRotating
export isinertial, isrotating, isbodycentric, isbarycentric
export @frame

# Transforms
export AstrodynamicsTransform
export Transform


using Unitful
using LabelledArrays
using DocStringExtensions
using CoordinateTransformations

include("Frames.jl")
include("Transforms.jl")

end # module
