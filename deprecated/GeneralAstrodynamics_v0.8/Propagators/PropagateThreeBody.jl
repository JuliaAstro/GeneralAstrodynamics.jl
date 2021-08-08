#
#   PropagateThreeBody.jl
#
#   Includes functions and structures for propagating orbits 
#   within the circular restricted three-body problem.
#





"""
Uses `DifferentialEquations` solvers to propagate normalized, 
Synodic CR3BP states Δt into the future. All keyword arguments are passed 
directly to `DifferentialEquations` solvers.
"""
propagate(orbit::NormalizedSynodicCR3BPOrbit, Δt::Unitful.Time; kwargs...) = propagate(
    orbit, nondimensionalize(Δt, normalized_time_unit(orbit.system)); kwargs...
)

