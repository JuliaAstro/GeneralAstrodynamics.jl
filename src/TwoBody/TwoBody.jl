""" 
    TwoBody

Provides structures & functions for the two-body problem.
"""
module TwoBody

# Dependencies 

using ..AbstractOrbits

using Reexport

using Logging
using Base: isapprox, isequal
using LinearAlgebra: ×, ⋅, norm
using DifferentialEquations
using ComponentArrays

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles

# Export data structures & constructors
export CartesianState, KeplerianState, TwoBodyState,
       TwoBodyOrbit, AbstractConic, Circular, 
       Elliptical, Parabolic, Hyperbolic

# Export functions
export  semimajor_axis, eccentricity, eccentricity_vector, inclination, true_anomoly, 
        periapsis_radius, apoapsis_radius, periapsis_velocity, apoapsis_velocity,      
        instantaneous_radius, instantaneous_velocity, orbital_period, time_since_periapsis, 
        mean_motion, mean_motion_vector, semi_parameter, eccentric_anomoly,
        specific_angular_momentum_vector, specific_angular_momentum, specific_energy,  
        isapprox, isequal, propagate, PropagationResult, kepler, TwoBodyOrbit, conic,
        KeplerianElements

# Include all module source code
include("TwoBodyStates.jl")
include("TwoBodyCalculations.jl")
include("TwoBodyPropagator.jl")
include("Kepler.jl")

end 