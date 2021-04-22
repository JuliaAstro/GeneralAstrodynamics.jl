"""
Contains structures, functions, and problem definitions
for core astrodynamics problems, including the 
Restricted Two-body Problem, the Circular Restricted
Three-body Problem, and the N-body problem.
"""
module Orbits

# Common stuctures and functions
export AbstractUnitfulStructure, AbstractState, AbstractFrame, AbstractSystem, AbstractOrbit, CartesianState
export MassParameter, lengthunit, timeunit, velocityunit, massparameterunit, coordinateframe
export position_vector, velocity_vector, scalar_position, scalar_velocity, epoch
export AbstractFrame, Inertial, Synodic, Perifocal, convert, epoch
export NormalizedLengthUnit, NormalizedTimeUnit, Trajectory, Manifold
export convert, show, eltype, isapprox, isequal

# Core data structures and functions for R2BP calculations 
export AbstractConic, Circular, Elliptical, Parabolic, Hyperbolic, Invalid
export KeplerianState, RestrictedTwoBodySystem, RestrictedTwoBodyOrbit
export R2BPState, R2BPSystem, R2BPOrbit, Orbit, CartesianOrbit, KeplerianOrbit
export mass_parameter, primary_synodic_position, secondary_synodic_position
export eccentricity, semimajor_axis, inclination, RAAN, argument_of_periapsis, true_anomoly
export conic, keplerian, cartesian, perifocal, semi_parameter, periapsis_radius
export specific_angular_momentum_vector, specific_angular_momentum, eccentricity_vector
export specific_energy, specific_potential_energy, mean_motion, mean_motion_vector
export eccentric_anomoly, time_since_periapsis, period
export kepler, lambert
export Sun, Mercury, Venus, Earth, Moon, Luna, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto

# Core data structures and functions for CR3BP calculations 
export CR3BPFrames, CR3BPSystem, CR3BPOrbit, CR3BPState
export SynodicCR3BPOrbit, NormalizedSynodicCR3BPOrbit, NormalizedSynodicSTMCR3BPOrbit
export SynodicCartesianState, InertialCartesianState, SynodicCartesianSTMState
export NormalizedCartesianState, MinimalCircularRestrictedThreeBodySystem
export CircularRestrictedThreeBodySystem, CircularRestrictedThreeBodyOrbit
export normalized_length_unit, normalized_time_unit, normalized_mass_parameter
export mass_parameters, primary_mass_parameter, secondary_mass_parameter
export time_scale_factor, nondimensionalize, redimensionalize
export nondimensionalize_length, nondimensionalize_time, nondimensionalize_velocity
export redimensionalize_length, redimensionalize_time, redimensionalize_velocity
export normalize, lagrange, inertial, synodic, accel, accel!, analyticalhalo
export potential_energy, jacobi_constant, zerovelocity_curves
export closest_approach, optimal_approach
export transform, transform_to_primary, transform_to_secondary
export SunEarth, EarthMoon

# Module Dependencies
using Reexport
using Contour
using StaticArrays
using LinearAlgebra
using Roots: find_zero

@reexport using Unitful, UnitfulAstro, UnitfulAngles

# Source Code 
include("Common/CommonTypes.jl")

include("R2BP/R2BPStates.jl")
include("R2BP/R2BPCalculations.jl")
include("R2BP/Kepler.jl")
include("R2BP/Lambert.jl")
include("R2BP/R2BPSystems.jl")

include("CR3BP/CR3BPStates.jl")
include("CR3BP/CR3BPCalculations.jl")
include("CR3BP/CR3BPSystems.jl")
include("CR3BP/CR3BPTransfers.jl")

end
