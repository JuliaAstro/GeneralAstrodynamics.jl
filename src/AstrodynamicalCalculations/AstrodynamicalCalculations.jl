"""
A module which provides common astrodynamics
calculations.

# Extended help

**Exports**

$(EXPORTS)

**Imports**

$(IMPORTS)
"""
module AstrodynamicalCalculations

export Circular, Elliptical, Parabolic, Hyperbolic
export conic, keplerian, cartesian
export perifocal, semimajor_axis
export specific_angular_momentum_vector, specific_angular_momentum
export specific_energy, C3, v_infinity, specific_potential_energy
export eccentricity_vector, eccentricity, semi_parameter
export distance, speed, periapsis_radius, apoapsis_radius, period
export true_anomoly, mean_motion, time_since_periapsis
export hohmann, SOI, SOA

using DocStringExtensions
using Unitful

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)

    $(DOCSTRING)
    """

@template (TYPES, CONSTANTS) =
    """
    $(TYPEDEF)

    $(DOCSTRING)
    """

include("R2BP/Conics.jl")
include("R2BP/R2BPCalculations.jl")

using ..AstrodynamicalStates
include("Hooks/AstrodynamicalStates.jl")

end # module
