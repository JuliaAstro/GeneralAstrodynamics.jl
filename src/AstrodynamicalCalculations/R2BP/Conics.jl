#
# Definitions for conic sections within 
# R2BP dynamics.
#

"""
$(TYPEDEF)

An abstract type for all conic sections.
"""
abstract type ConicSection end

"""
$(TYPEDEF)

An abstract type representing the Circular 
conic section.
"""
abstract type Circular <: ConicSection end

"""
$(TYPEDEF)

An abstract type representing the Elliptical 
conic section.
"""
abstract type Elliptical <: ConicSection end

"""
$(TYPEDEF)

An abstract type representing the Parabolic 
conic section.
"""
abstract type Parabolic <: ConicSection end

"""
$(TYPEDEF)

An abstract type representing the Hyperbolic 
conic section.
"""
abstract type Hyperbolic <: ConicSection end
