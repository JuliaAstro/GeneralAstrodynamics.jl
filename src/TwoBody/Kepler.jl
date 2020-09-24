#
#   Kepler.jl
#
#   Solves Kepler's problem for TwoBody orbits.
#

function kepler(orbit::AbstractOrbit, Δt=orbital_period(orbit))

# Find select orbital elements
a = semimajor_axis(orbit)
e = eccentricity(orbit)

# Guess χₙ
# if iselliptical(orbit) && isrevolution(orbit, Δt)
    

# end

end

"""
    isrevolution(orbit::TwoBodyOrbit, Δt::Unitful.Time)

Checks if Δt is a multiple of `orbit`'s orbital period.
"""
isrevolution(orbit::TwoBodyOrbit, Δt::Unitful.Time) = mod(Δt, orbital_period(orbit)) ≈ orbital_period(orbit)