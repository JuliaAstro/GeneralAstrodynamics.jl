#
# Maneuver calculations for orbits within the Twobody Problem
#

"""
    escape_radius(r₀, v₀, aₜ)
    escape_radius(orbit::Orbit, m::ConstantManeuver)

Provides the radius of escape.
"""
escape_radius(r₀, v₀, aₜ) = r₀ * v₀ / (20 * aₜ^2 * r₀^2)^(1/4)
escape_radius(orbit::Orbit, m::ConstantManeuver) = escape_radius(radius(orbit), velocity(orbit), m.aₜ)

"""
    escape_velocity(r₀, v₀, aₜ, μ)
    escape_velocity(orbit::Orbit, m::ConstantManeuver)

Provides the velocity at escape.
"""
escape_velocity(r₀, v₀, aₜ, μ) = √(2 * μ / escape_radius(r₀, v₀, aₜ))
escape_velocity(orbit::Orbit, m::ConstantManeuver) = escape_velocity(radius(orbit), velocity(orbit), m.aₜ, orbit.body.μ)