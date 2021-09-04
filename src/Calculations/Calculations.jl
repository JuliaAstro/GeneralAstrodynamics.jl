"""
A module which provides common astrodynamics
calculations.

# Extended help

**Exports**

$(EXPORTS)

**Imports**

$(IMPORTS)
"""
module Calculations

export Circular, Elliptical, Parabolic, Hyperbolic
export conic, keplerian, cartesian
export perifocal, semimajor_axis
export kepler, lambert, lambert_universal, lambert_lancaster_blanchard
export specific_angular_momentum_vector, specific_angular_momentum
export specific_energy, C3, v_infinity, specific_potential_energy
export eccentricity_vector, eccentricity, semi_parameter
export distance, speed, periapsis_radius, apoapsis_radius, period
export true_anomoly, mean_motion, time_since_periapsis
export hohmann, SOI, SOA, normalize, redimension, analyticalhalo
export jacobi_constant, zerovelocity_curves

using DocStringExtensions
using Unitful
using Requires
using Contour
using StaticArrays
using AstroTime
using LinearAlgebra

import Dates
import Roots: find_zero

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

Unitful.@derived_dimension MassParameter Unitful.𝐋^3/Unitful.𝐓^2

@doc """
A unit dimension alias for length^3 / time^2. This is a common
dimension used in astrodynamics calculations.
"""
MassParameter
    
include(joinpath("R2BP",  "Conics.jl"))
include(joinpath("R2BP",  "R2BPCalculations.jl"))
include(joinpath("R2BP",  "Kepler.jl"))
include(joinpath("R2BP",  "Lambert.jl"))
include(joinpath("CR3BP", "CR3BPCalculations.jl"))

end # module
