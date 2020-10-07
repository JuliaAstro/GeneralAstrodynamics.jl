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

# Newton's Gravitation Constant
import PhysicalConstants.CODATA2018
G = 1.0 * CODATA2018.G

@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro, UnitfulAngles
@reexport using Plots, Plots.PlotMeasures

# Export data structures, constants, and constructors
export TwoBodySystem, TwoBodyOrbit, AbstractConic, Circular, 
       Elliptical, Parabolic, Hyperbolic, Body, G, CelestialBody,
       Sun, Mercury, Venus, Earth, Moon, Luna, Mars, Jupiter, 
       Saturn, Uranus, Neptune, Pluto

# Export functions
export  semimajor_axis, eccentricity, eccentricity_vector, inclination, true_anomoly, 
        periapsis_radius, apoapsis_radius, periapsis_velocity, apoapsis_velocity,      
        radius, radius_vector, velocity, velocity_vector, orbital_period, 
        time_since_periapsis, mean_motion, mean_motion_vector, semi_parameter, conic_anomoly,
        specific_angular_momentum_vector, specific_angular_momentum, specific_energy,  
        isapprox, isequal, propagate_twobody, TwobodyPropagationResult, kepler, conic,
        orbital_elements, cartesian, isinvalid, twobody_plot3d

# Include all module source code
include("twobody_states.jl")
include("twobody_calculations.jl")
include("propagate_twobody.jl")
include("kepler.jl")
include("plot_twobody.jl")

end 