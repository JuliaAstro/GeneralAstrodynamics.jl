#
# Maneuver calculations for orbits within the Twobody Problem
#

"""
Provides the radius of escape.
"""
escape_radius(r₀, v₀, aₜ) = r₀ * v₀ / (20 * aₜ^2 * r₀^2)^(1/4)
escape_radius(orbit::T, m::ConstantManeuver) where T <: RestrictedTwoBodySystem = escape_radius(radius(orbit), velocity(orbit), m.aₜ)

"""
Provides the velocity at escape.
"""
escape_velocity(r₀, v₀, aₜ, μ) = √(2 * μ / escape_radius(r₀, v₀, aₜ))
escape_velocity(orbit::T, m::ConstantManeuver) where T <: RestrictedTwoBodySystem = escape_velocity(radius(orbit), velocity(orbit), m.aₜ, orbit.body.μ)

"""
Provides time delta from the provided initial orbit to escape.
"""
escape_time(r₀, v₀, aₜ) = upreferred(v₀ / aₜ) * (1 - (upreferred(20aₜ^2 * r₀^2) / upreferred(v₀^4))^(1/8))
escape_time(orbit::T, m::ConstantManeuver) where T <: RestrictedTwoBodySystem = escape_time(radius(orbit), velocity(orbit), m.aₜ)

"""
Provides the path length from the initial condition to escape.
"""
escape_path_length(r₀, v₀, aₜ) = upreferred(v₀^2 / 2aₜ) * (1 - upreferred((1/v₀) * (20aₜ^2 * r₀^2)^(1/4)))
escape_path_length(orbit::T, m::ConstantManeuver) where T <: RestrictedTwoBodySystem = escape_path_length(radius(orbit), velocity(orbit), m.aₜ)