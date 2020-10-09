""" 
    TwoBody

Provides structures & functions for the two-body problem.
"""
module TwoBody

# Dependencies 

using ..AbstractTypes

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using DifferentialEquations
using ComponentArrays
using StaticArrays

# Newton's Gravitation Constant
import PhysicalConstants.CODATA2018
G = 1.0 * CODATA2018.G

@reexport using Unitful, UnitfulAstro, UnitfulAngles

# Export data structures, constants, and constructors
export TwoBodySystem, Orbit, AbstractConic, Circular, InvalidOrbit,
       Elliptical, Parabolic, Hyperbolic, Body, CelestialBody,
       Sun, Mercury, Venus, Earth, Moon, Luna, Mars, Jupiter, 
       Saturn, Uranus, Neptune, Pluto

# Export functions
export  semimajor_axis, eccentricity, eccentricity_vector, inclination, true_anomoly, 
        periapsis_radius, apoapsis_radius, periapsis_velocity, apoapsis_velocity,      
        radius, radius_vector, velocity, velocity_vector, orbital_period, 
        time_since_periapsis, mean_motion, mean_motion_vector, semi_parameter, conic_anomoly,
        specific_angular_momentum_vector, specific_angular_momentum, specific_energy,  
        isapprox, isequal, TwobodyPropagationResult, kepler, conic,
        orbital_elements, cartesian, isinvalid

# Include all module source code
include("twobody_states.jl")
include("twobody_calculations.jl")
include("kepler.jl")

end 