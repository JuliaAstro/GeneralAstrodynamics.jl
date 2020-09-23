""" 
    TwoBody

Provides structures & functions for the two-body problem.
"""
module TwoBody

# Dependencies 

using ..AbstractOrbits

using Reexport

using Base: isapprox, isequal
using Logging
using LinearAlgebra: ×, ⋅, norm
using DifferentialEquations
using ComponentArrays

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles

# Export data structures & constructors
export Body, earth, sun
export CartesianState, KeplerianState

# Export functions
export  semimajor_axis, eccentricity, eccentricity_vector, inclination, true_anomoly, 
        periapsis_radius, apoapsis_radius, periapsis_velocity, apoapsis_velocity,      
        instantaneous_radius, instantaneous_velocity, orbital_period, time_since_periapsis, 
        mean_motion, mean_motion_vector, semi_parameter, eccentric_anomoly,
        specific_angular_momentum_vector, specific_angular_momentum, specific_energy,  
        isapprox, isequal, propagate, PropagationResult

# Include all module source code
include("TwoBodyStates.jl")
include("TwoBodyCalculations.jl")
include("TwoBodyPropagator.jl")

end 