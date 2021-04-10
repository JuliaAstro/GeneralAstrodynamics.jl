""" 
Provides structures & functions for the two-body problem.
"""
module TwoBody

# Dependencies 

using Reexport
using ..AstrodynamicsCore

using Unitful, UnitfulAngles

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

using Crayons
using LinearAlgebra: norm, cross, ×, dot, ⋅
using StaticArrays

# Newton's Gravitation Constant
import PhysicalConstants.CODATA2018
G = 1.0 * CODATA2018.G

# Export data structures, constants, and constructors
export RestrictedTwoBodyOrbit, RestrictedTwoBodySystem, KeplerianState, AbstractConic, Circular,
       Elliptical, Parabolic, Hyperbolic, Invalid, Body, CelestialBody,
       Sun, Mercury, Venus, Earth, Moon, Luna, Mars, Jupiter, 
       Saturn, Uranus, Neptune, Pluto, G

# Export functions
export  semimajor_axis, semi_parameter, eccentricity, 
        eccentricity_vector, inclination, true_anomoly, 
        periapsis_radius, apoapsis_radius, periapsis_velocity, 
        apoapsis_velocity, scalar_position, scalar_velocity, period, 
        position_vector, velocity_vector, Orbit, perifocal_radius,
        mass, mass_parameter, perifocal, RAAN, argument_of_periapsis,
        time_since_periapsis, mean_motion, mean_motion_vector, 
        eccentric_anomoly, specific_angular_momentum_vector, 
        specific_angular_momentum, specific_energy, specific_potential_energy,
        isapprox, isequal, TwobodyPropagationResult, kepler, lambert,
        conic, keplerian, cartesian, promote, convert, Float16, Float32, Float64,
        BigFloat, CartesianOrbit, KeplerianOrbit, R2BP

# Include all module source code
include("TwoBodyStates.jl")
include("TwoBodyCalculations.jl")
include("Kepler.jl")
include("Lambert.jl")

include("TwoBodySystems.jl")

@reexport using .Systems

end 