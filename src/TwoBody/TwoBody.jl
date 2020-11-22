""" 
Provides structures & functions for the two-body problem.
"""
module TwoBody

include("../Misc/DocStringExtensions.jl")

# Dependencies 

using Reexport
@reexport using ..CommonTypes

using Crayons
using LinearAlgebra: norm, cross, ×, dot, ⋅
using StaticArrays: SVector, @SVector, SMatrix, @SMatrix

# Newton's Gravitation Constant
import PhysicalConstants.CODATA2018
G = 1.0 * CODATA2018.G

# Export data structures, constants, and constructors
export TwoBodySystem, Orbit, AbstractConic, Circular, InvalidOrbit,
       Elliptical, Parabolic, Hyperbolic, Invalid, Body, CelestialBody,
       Sun, Mercury, Venus, Earth, Moon, Luna, Mars, Jupiter, 
       Saturn, Uranus, Neptune, Pluto, G

# Export functions
export  semimajor_axis, semi_parameter, eccentricity, 
        eccentricity_vector, inclination, true_anomoly, 
        periapsis_radius, apoapsis_radius, periapsis_velocity, 
        apoapsis_velocity, radius, velocity, orbital_period, 
        mass, mass_parameter, perifocal, inertial,
        time_since_periapsis, mean_motion, mean_motion_vector, 
        conic_anomoly, specific_angular_momentum_vector, 
        specific_angular_momentum, specific_energy,  
        isapprox, isequal, TwobodyPropagationResult, kepler, lambert,
        conic, keplerian, cartesian, isinvalid, promote, convert

# Include all module source code
include("twobody_states.jl")
include("twobody_calculations.jl")
include("kepler.jl")
include("lambert.jl")

end 