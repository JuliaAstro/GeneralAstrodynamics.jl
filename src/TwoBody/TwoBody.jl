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

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles

# Export data structures & constructors
export TwoBodySystem, TwoBodyOrbit, AbstractConic, Circular, 
       Elliptical, Parabolic, Hyperbolic, Body, Earth, Sun

# Export functions
export  semimajor_axis, eccentricity, eccentricity_vector, inclination, true_anomoly, 
        periapsis_radius, apoapsis_radius, periapsis_velocity, apoapsis_velocity,      
        radius, radius_vector, velocity, velocity_vector, orbital_period, 
        time_since_periapsis, mean_motion, mean_motion_vector, semi_parameter, conic_anomoly,
        specific_angular_momentum_vector, specific_angular_momentum, specific_energy,  
        isapprox, isequal, propagate_twobody, TwobodyPropagationResult, kepler, conic,
        orbital_elements, cartesian

# Include all module source code
include("states.jl")
include("calculations.jl")
include("propagator.jl")
include("kepler.jl")

end 