#
# Data structures and types for orbit maneuvers
#

"""
    AbstractManeuver

An abstract type for all orbit maneuvers.
"""
abstract type AbstractManeuver end

"""
    TwoBodyManeuver

An abstract type for all twobody maneuvers.
"""
abstract type TwoBodyManeuver<:AbstractManeuver end

"""
    ConstantManeuver

A type for constant, continuous thrust.
"""
struct ConstantManeuver<:TwoBodyManeuver 

    aₜ::Unitful.Acceleration

end