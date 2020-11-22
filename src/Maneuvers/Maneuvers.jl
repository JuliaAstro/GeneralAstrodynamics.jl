"""
Provides calculations for orbit maneuvers.
"""
module Maneuvers

include("../Misc/DocStringExtensions.jl")

using Reexport
@reexport using ..CommonTypes

using ..NBody
using ..TwoBody
using LinearAlgebra: norm, cross, ×, dot, ⋅

export AbstractManeuver, TwoBodyManeuver, ConstantManeuver
export escape_radius, escape_velocity, escape_time, escape_path_length

include("maneuver_types.jl")
include("twobody_maneuver_calculations.jl")

end