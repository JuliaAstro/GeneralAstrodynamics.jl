#
# Data structures and types for orbit maneuvers
#

"""
    AbstractManeuver

An abstract type for all orbit maneuvers.
"""
abstract type AbstractManeuver end

"""
    ConstantManeuver

A type for constant, continuous thrust.
"""
struct ConstantManeuver<:AbstractManeuver 

    aₜ::Unitful.Acceleration

end