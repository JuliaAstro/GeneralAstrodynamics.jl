"""
    Maneuvers

Provides calculations for orbit maneuvers.
"""
module Maneuvers

using ..CommonTypes
using ..NBody
using ..TwoBody

using Reexport

using Logging

export AbstractManeuver, ConstantManeuver
export escape_radius, escape_velocity, escape_time

include("maneuver_types.jl")
include("twobody_maneuver_calculations.jl")

end