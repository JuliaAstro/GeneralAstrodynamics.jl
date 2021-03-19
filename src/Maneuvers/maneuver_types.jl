#
# Data structures and types for orbit maneuvers
#

"""
An abstract type for all orbit maneuvers.
"""
abstract type AbstractManeuver end

"""
An abstract type for all twobody maneuvers.
"""
abstract type TwoBodyManeuver <: AbstractManeuver end

"""
A type for constant, continuous thrust.
"""
struct ConstantManeuver <: TwoBodyManeuver 

    aₜ::Unitful.Acceleration

end