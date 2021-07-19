#
# Provides calculations directory on `Orbit` types.
#

#
# Defines AstrodynamicalCalculations methods for 
# AstrodynamicalStates types
#

"""
Returns the position vector of the `CartesianState`.
"""
Base.position(state::CartesianState) = AstrodynamicalStates.get_r(state) * lengthunit(state)

"""
Returns the scalar position of the `CartesianState`.
"""
distance(state::CartesianState) = norm(position(state))

"""
Returns the velocity vector for the `CartesianState`.
"""
velocity(state::CartesianState) = AstrodynamicalStates.get_v(state) * velocityunit(state)

"""
Returns the scalar velocity of the `CartesianState`.
"""
speed(state::CartesianState) = norm(velocity(state))

"""
Returns orbital eccentricity, `e`.
"""
eccentricity(state::KeplerianState) = AstrodynamicalStates.get_e(state)

"""
Returns semi-major axis, `a`.
"""
semimajor_axis(state::KeplerianState) = AstrodynamicalStates.get_a(state) * lengthunit(state)

"""
Returns orbital inclination, `i`.
"""
inclination(state::KeplerianState) = AstrodynamicalStates.get_i(state) * angularunit(state)

"""
Returns the right ascension of the ascending node (R.A.A.N), `Ω`.
"""
RAAN(state::KeplerianState) = AstrodynamicalStates.get_Ω(state) * angularunit(state)

"""
Returns the argument of periapsis, `ω`.
"""
argument_of_periapsis(state::KeplerianState) = AstrodynamicalStates.get_ω(state) * angularunit(state)

"""
Returns the true anomoly, `ν`.
"""
true_anomoly(state::KeplerianState) = AstrodynamicalStates.get_ν(state) * angularunit(state)

#=

# This isn't necessary yet, since they're all in the same package!

conic(orbit::T) where T<:RestrictedTwoBodyOrbit = conic(eccentricity(orbit))

keplerian(rᵢ, vᵢ, body::RestrictedTwoBodySystem) = keplerian(rᵢ, vᵢ, mass_parameter(body))
keplerian(orbit::KeplerianOrbit) = eccentricity(orbit.state), semimajor_axis(orbit.state), inclination(orbit.state), RAAN(orbit.state), argument_of_periapsis(orbit.state), true_anomoly(orbit.state)
keplerian(orbit::CartesianOrbit) = keplerian(position_vector(orbit.state), velocity_vector(orbit.state), orbit.system)

cartesian(e, a, i, Ω, ω, ν, body::RestrictedTwoBodySystem) = cartesian(e, a, i, Ω, ω, ν, mass_parameter(body))
cartesian(orbit::CartesianOrbit) = position_vector(orbit.state), velocity_vector(orbit.state)
cartesian(orbit::KeplerianOrbit) = cartesian(eccentricity(orbit.state), semimajor_axis(orbit.state), inclination(orbit.state), RAAN(orbit.state), argument_of_periapsis(orbit.state), true_anomoly(orbit.state), orbit.body)

function perifocal(orbit::T) where T <: RestrictedTwoBodyOrbit
    return perifocal(
        inclination(orbit),
        RAAN(orbit),
        argument_of_periapsis(orbit),
        position_vector(orbit),
        velocity_vector(orbit)
    )
end

semimajor_axis(orbit::CartesianOrbit) = semimajor_axis(scalar_position(orbit), scalar_velocity(orbit), mass_parameter(orbit.system)) 
semimajor_axis(orbit::KeplerianOrbit) = orbit.state.a * lengthunit(orbit.state)

specific_angular_momentum_vector(orbit::CartesianOrbit) = specific_angular_momentum_vector(position_vector(orbit.state), velocity_vector(orbit.state))
specific_angular_momentum_vector(orbit::KeplerianOrbit) = specific_angular_momentum_vector(CartesianOrbit(orbit))

specific_angular_momentum(orbit::CartesianOrbit) = specific_angular_momentum(position_vector(orbit.state), velocity_vector(orbit.state))
specific_angular_momentum(orbit::KeplerianOrbit) = specific_angular_momentum(cartesian(orbit)...)

specific_energy(orbit::CartesianOrbit) = specific_energy(scalar_position(orbit.state), scalar_velocity(orbit.state), mass_parameter(orbit.system))
specific_energy(orbit::KeplerianOrbit) = specific_energy(semimajor_axis(orbit.state), mass_parameter(orbit.system))

C3(orbit::RestrictedTwoBodyOrbit) = C3(scalar_position(orbit), scalar_velocity(orbit), mass_parameter(orbit.system))

v_infinity(orbit::RestrictedTwoBodyOrbit) = v_infinity(scalar_position(orbit), scalar_velocity(orbit), mass_parameter(orbit.system))

eccentricity_vector(orbit::CartesianOrbit) = eccentricity_vector(position_vector(orbit.state), velocity_vector(orbit.state), mass_parameter(orbit.system))
eccentricity_vector(orbit::KeplerianOrbit) = eccentricity_vector(CartesianOrbit(orbit))

eccentricity(orbit::CartesianOrbit) = eccentricity(position_vector(orbit.state), velocity_vector(orbit.state), mass_parameter(orbit.system))
eccentricity(orbit::KeplerianOrbit) = orbit.state.e

semi_parameter(orbit::RestrictedTwoBodyOrbit) = semi_parameter(semimajor_axis(orbit), eccentricity(orbit))

scalar_position(orbit::KeplerianOrbit) = scalar_position(semi_parameter(orbit), eccentricity(orbit), true_anomoly(orbit))
scalar_position(orbit::CartesianOrbit) = norm(position_vector(orbit))

scalar_position(orbit::KeplerianOrbit) = scalar_position(semi_parameter(orbit), eccentricity(orbit), true_anomoly(orbit))
scalar_position(orbit::CartesianOrbit) = norm(position_vector(orbit))

position_vector(orbit::CartesianOrbit) = position_vector(orbit.state)
position_vector(orbit::KeplerianOrbit) = position_vector(CartesianOrbit(orbit))

scalar_velocity(orbit::KeplerianOrbit) = scalar_velocity(scalar_position(orbit), semimajor_axis(orbit), mass_parameter(orbit.system))
scalar_velocity(orbit::CartesianOrbit) = norm(velocity_vector(orbit))

velocity_vector(orbit::CartesianOrbit) = velocity_vector(orbit.state)
velocity_vector(orbit::KeplerianOrbit) = velocity_vector(CartesianOrbit(orbit))

periapsis_radius(orbit::T) where T<:RestrictedTwoBodyOrbit = periapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

apoapsis_radius(orbit::T) where T<:RestrictedTwoBodyOrbit = apoapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

periapsis_velocity(orbit::T) where T<:RestrictedTwoBodyOrbit = scalar_velocity(periapsis_radius(orbit), semimajor_axis(orbit), mass_parameter(orbit.system))

apoapsis_velocity(orbit::T) where T<:RestrictedTwoBodyOrbit = scalar_velocity(apoapsis_radius(orbit), semimajor_axis(orbit), mass_parameter(orbit.system))

"""
Returns mass `m`.
"""
mass(body::RestrictedTwoBodySystem) = mass_parameter(body) / G

"""
Returns mass parameter `μ`.
"""
mass_parameter(body::RestrictedTwoBodySystem) = body.μ * massparameterunit(body)

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{Parabolic}) = Inf * timeunit(orbit.state)

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{Hyperbolic}) = NaN * timeunit(orbit.state)

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{Invalid}) = NaN * timeunit(orbit.state)

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{C}) where C<:Union{Elliptical, Circular} = period(semimajor_axis(orbit), mass_parameter(orbit.system))

true_anomoly(orbit::KeplerianOrbit) = orbit.state.ν * angularunit(orbit.state)

true_anomoly(orbit::CartesianOrbit) = true_anomoly(scalar_position(orbit), specific_angular_momentum(orbit), eccentricity(orbit), mass_parameter(orbit.system))

mean_motion(orbit::T) where T<:RestrictedTwoBodyOrbit = mean_motion(semimajor_axis(orbit), mass_parameter(orbit.system))

"""
Returns mean motion vector, n̄.
"""
function mean_motion_vector(orbit::T) where T<:RestrictedTwoBodyOrbit
#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}(0, 0, 1)
    return k̂ × specific_angular_momentum_vector(orbit)
end


"""
Returns eccentric anomoly, E, parabolic anomoly, B, or hyperbolic 
anomoly, H. 
"""
function eccentric_anomoly(orbit::T) where T <: RestrictedTwoBodyOrbit
    e = eccentricity(orbit)
    ν = true_anomoly(orbit)
    return acos(u"rad", (e + cos(ν) / (1 + e * cos(ν))))
end

"""
Returns time since periapsis, t.
"""
time_since_periapsis(orbit::T) where T <: RestrictedTwoBodyOrbit = time_since_periapsis(mean_motion(orbit), eccentricity(orbit), eccentric_anomoly(orbit))

"""
Returns orbital inclination, i.
"""
inclination(orbit::KeplerianOrbit) =  orbit.state.i * angularunit(orbit.state)

"""
Returns orbital inclination, i.
"""
inclination(orbit::CartesianOrbit) = inclination(KeplerianOrbit(orbit))

"""
Returns Right Ascension of the Ascending Node, Ω.
"""
RAAN(orbit::KeplerianOrbit) = orbit.state.Ω * angularunit(orbit.state)

"""
Returns Right Ascension of the Ascending Node, Ω.
"""
RAAN(orbit::CartesianOrbit) = RAAN(KeplerianOrbit(orbit))

"""
Returns the Argument of Periapsis, ω.
"""
argument_of_periapsis(orbit::KeplerianOrbit) = orbit.state.ω * angularunit(orbit.state)

"""
Returns the Argument of Periapsis, ω.
"""
argument_of_periapsis(orbit::CartesianOrbit) = argument_of_periapsis(KeplerianOrbit(orbit))


"""
Returns true if all elements in each system are within `atol` of the other.
"""
function Base.isapprox(c1::RestrictedTwoBodyOrbit, c2::RestrictedTwoBodyOrbit; atol=1e-6)

    return all(ustrip.(abs.(position_vector(c1) - position_vector(c2))) .< atol) &&
           all(ustrip.(abs.(velocity_vector(c1) - velocity_vector(c2))) .< atol) &&
           ustrip(upreferred(abs(eccentricity(c1) - eccentricity(c2)))) < atol &&
           ustrip(upreferred(abs(semimajor_axis(c1) - semimajor_axis(c2)))) < atol &&
           ustrip(upreferred(abs(mod(inclination(c1), 180u"°") - mod(inclination(c2), 180u"°")))) < atol &&
           ustrip(upreferred(abs(mod(RAAN(c1), 360u"°") - mod(RAAN(c2), 360u"°")))) < atol &&
           ustrip(upreferred(abs(mod(argument_of_periapsis(c1), 360u"°") - mod(argument_of_periapsis(c2), 360u"°")))) < atol &&
           ustrip(upreferred(abs(mod(true_anomoly(c1), 360u"°") - mod(true_anomoly(c2), 360u"°")))) < atol &&
           isapprox(c1.system, c2.system; atol=atol)

end

"""
Returns true if all elements of each system are identically equal.
"""
function Base.isequal(c1::RestrictedTwoBodyOrbit, c2::RestrictedTwoBodyOrbit)

    return all(position_vector(c1) .== position_vector(c2)) &&
           all(velocity_vector(c1) .== velocity_vector(c2)) &&
           eccentricity(c1) == eccentricity(c2) &&
           semimajor_axis(c1) == semimajor_axis(c2) &&
           mod(inclination(c1), 180u"°") == mod(inclination(c2), 180u"°") &&
           mod(RAAN(c1), 360u"°") == mod(RAAN(c2), 360u"°") &&
           mod(argument_of_periapsis(c1), 360u"°") == mod(argument_of_periapsis(c2), 360u"°") &&
           mod(true_anomoly(c1), 360u"°") == mod(true_anomoly(c2), 360u"°") &&
           c1.system == c2.system
           
end

specific_potential_energy(orbit::CartesianOrbit) = specific_potential_energy(scalar_position(orbit.state), mass_parameter(orbit.system))
specific_potential_energy(orbit::KeplerianOrbit) = specific_potential_energy(CartesianOrbit(orbit))

"""
Returns true if all elements are within `atol` of the other.
"""
function Base.isapprox(b1::RestrictedTwoBodySystem, b2::RestrictedTwoBodySystem; atol=1e-6)
    return ustrip(upreferred(abs(mass_parameter(b1) - mass_parameter(b2)))) < atol
end

"""
Returns true if all elements are identically equal.
"""
function Base.isequal(b1::RestrictedTwoBodySystem, b2::RestrictedTwoBodySystem)
    return mass_parameter(b1) == mass_parameter(b2)
end


=#