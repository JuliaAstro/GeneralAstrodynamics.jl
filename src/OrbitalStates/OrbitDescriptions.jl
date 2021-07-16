#
# Orbit descriptions.
#

"""
$(TYPEDEF)

A supertype for all astrodynamical systems.
"""
abstract type AstrodynamicalSystem{F, LU, TU, AU, P} end

#=
    1) Make a `UnitfulLabelledArrays` type
    2) 
=#

"""
$(TYPEDEF)

A Restricted Two-body Problem system.
"""


"""
$(TYPEDEF)

A supertype for all single-point orbit descriptions.
"""
abstract type OrbitalState{F, LU, TU, AU, S} end

"""
$(TYPEDEF) 

A supertype for all bodies in space.
"""
abstract type OrbitBody end

