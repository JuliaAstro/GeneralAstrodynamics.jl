#
# Provides calculations directory on `Orbit` types.
#

#
# Defines AstrodynamicalCalculations methods for 
# AstrodynamicalStates types
#

using Lazy 

"""
Converts a `CartesianState` to a `KeplerianState`.
"""
keplerian(state::CartesianState, μ) = keplerian(position(state), velocity(state), μ)

"""
Converts a `CartesianR2BPOrbit` to a `KeplerianR2BPOrbit`.
"""
keplerian(orbit::CartesianR2BPOrbit) = keplerian(state(orbit), massparameter(system(orbit)))

"""
No-operation applied, orbit is already `Keplerian`.
"""
keplerian(orbit::KeplerianR2BPOrbit) = orbit

"""
Converts a `KeplerianState` to a `CartesianState`.
"""
cartesian(state::KeplerianState, μ) = cartesian(
    eccentricity(state),
    semimajor_axis(state),
    inclination(state),
    RAAN(state),
    argument_of_periapsis(state),
    true_anomoly(state),
    μ
)

"""
Converts a `KeplerianR2BPOrbit` to a `CartesianR2BPOrbit`.
"""
cartesian(orbit::KeplerianR2BPOrbit) = cartesian(state(orbit), massparameter(system(orbit)))

"""
No-operation applied, orbit is already `Cartesian`.
"""
cartesian(orbit::CartesianR2BPOrbit) = orbit

"""
Returns the position vector of the `CartesianState`.
"""
Base.position(state::CartesianState) = AstrodynamicalStates.get_r(state) * lengthunit(state)

"""
Returns the position vector of the `CartesianR2BPOrbit`.
"""
Base.position(orbit::CartesianOrbit) = position(state(orbit))

"""
Returns the position vector of the `KeplerianR2BPOrbit`.
"""
Base.position(orbit::KeplerianR2BPOrbit) = position(cartesian(state(orbit), massparameter(state(orbit))))

"""
Returns the scalar position of the `CartesianState`.
"""
distance(state::CartesianState) = norm(position(state))

"""
Returns the scalar position of the `Orbit`.
"""
distance(orbit::Orbit) = norm(position(orbit))

"""
Returns the velocity vector for the `CartesianState`.
"""
velocity(state::CartesianState) = AstrodynamicalStates.get_v(state) * velocityunit(state)

"""
Returns the velocity vector for the `CartesianR2BPOrbit`.
"""
velocity(orbit::CartesianOrbit) = velocity(state(orbit))

"""
Returns the scalar velocity of the `CartesianState`.
"""
speed(state::CartesianState) = norm(velocity(state))

"""
Returns the scalar velocity of the `CartesianOrbit`.
"""
speed(orbit::CartesianOrbit) = norm(velocity(state(orbit)))

"""
Returns the scalar velocity of the `KeplerianState`.
"""
speed(orbit::KeplerianR2BPOrbit) = speed(semi_parameter(semimajor_axis(state(orbit)), eccentricity(state(orbit))), semimajor_axis(state(orbit)), massparameter(system(orbit)))

"""
Returns orbital eccentricity, `e`.
"""
eccentricity(state::KeplerianState) = AstrodynamicalStates.get_e(state)
eccentricity(orbit::KeplerianOrbit) = eccentricity(state(orbit))
eccentricity(orbit::CartesianOrbit) = eccentricity(position(state(orbit)), velocity(state(orbit)), massparameter(system(orbit)))

"""
Returns semi-major axis, `a`.
"""
semimajor_axis(state::KeplerianState) = AstrodynamicalStates.get_a(state) * lengthunit(state)
semimajor_axis(orbit::CartesianR2BPOrbit) = semimajor_axis(distance(state(orbit)), speed(state(orbit)), massparameter(system(orbit)))
semimajor_axis(orbit::KeplerianR2BPOrbit) = semimajor_axis(state(orbit))

"""
Returns orbital inclination, `i`.
"""
inclination(state::KeplerianState) = AstrodynamicalStates.get_i(state) * angularunit(state)
inclination(orbit::KeplerianR2BPOrbit) = inclination(state(orbit))

"""
Returns the right ascension of the ascending node (R.A.A.N), `Ω`.
"""
RAAN(state::KeplerianState) = AstrodynamicalStates.get_Ω(state) * angularunit(state)
RAAN(orbit::KeplerianR2BPOrbit) = RAAN(state(orbit))

"""
Returns the argument of periapsis, `ω`.
"""
argument_of_periapsis(state::KeplerianState) = AstrodynamicalStates.get_ω(state) * angularunit(state)
argument_of_periapsis(orbit::KeplerianR2BPOrbit) = argument_of_periapsis(state(orbit))

"""
Returns the true anomoly, `ν`.
"""
true_anomoly(state::KeplerianState) = AstrodynamicalStates.get_ν(state) * angularunit(state)
true_anomoly(orbit::KeplerianR2BPOrbit) = true_anomoly(state(orbit))

"""
Returns the conic section of the `Orbit`.
"""
conic(orbit::R2BPOrbit) = conic(eccentricity(orbit))

"""
Returns the specific angular momentum vector of the `CartesianR2BPOrbit`.
"""
specific_angular_momentum_vector(orbit::CartesianR2BPOrbit) = specific_angular_momentum_vector(position(state(orbit)), velocity_vector(state(orbit)))

"""
Returns the specific angular momentum of the `CartesianR2BPOrbit`.
"""
specific_angular_momentum(orbit::CartesianR2BPOrbit) = specific_angular_momentum(position(state(orbit)), velocity_vector(state(orbit)))

"""
Returns the specific energy.
"""
specific_energy(orbit::CartesianR2BPOrbit) = specific_energy(distance(state(orbit)), speed(state(orbit)), massparameter(system(orbit)))
specific_energy(orbit::KeplerianR2BPOrbit) = specific_energy(semimajor_axis(orbit), massparameter(system(orbit)))

"""
Returns the orbit's C3.
"""
C3(orbit::R2BPOrbit) = C3(distance(orbit), speed(orbit), massparameter(system(orbit)))

"""
Returns escape velocity, $v_\infty$.
"""
v_infinity(orbit::R2BPOrbit) = v_infinity(distance(orbit), speed(orbit), massparameter(system(orbit)))

eccentricity_vector(orbit::CartesianOrbit) = eccentricity_vector(position_vector(state(orbit)), velocity_vector(state(orbit)), massparameter(system(orbit)))
eccentricity_vector(orbit::KeplerianOrbit) = eccentricity_vector(CartesianOrbit(orbit))

eccentricity(orbit::CartesianOrbit) = eccentricity(position_vector(state(orbit)), velocity_vector(state(orbit)), massparameter(system(orbit)))
eccentricity(orbit::KeplerianOrbit) = state(orbit).e

semi_parameter(orbit::RestrictedTwoBodyOrbit) = semi_parameter(semimajor_axis(orbit), eccentricity(orbit))

distance(orbit::KeplerianOrbit) = distance(semi_parameter(orbit), eccentricity(orbit), true_anomoly(orbit))
distance(orbit::CartesianOrbit) = norm(position_vector(orbit))

distance(orbit::KeplerianOrbit) = distance(semi_parameter(orbit), eccentricity(orbit), true_anomoly(orbit))
distance(orbit::CartesianOrbit) = norm(position_vector(orbit))

position_vector(orbit::CartesianOrbit) = position_vector(state(orbit))
position_vector(orbit::KeplerianOrbit) = position_vector(CartesianOrbit(orbit))

speed(orbit::KeplerianOrbit) = speed(distance(orbit), semimajor_axis(orbit), massparameter(system(orbit)))
speed(orbit::CartesianOrbit) = norm(velocity_vector(orbit))

velocity_vector(orbit::CartesianOrbit) = velocity_vector(state(orbit))
velocity_vector(orbit::KeplerianOrbit) = velocity_vector(CartesianOrbit(orbit))

periapsis_radius(orbit::T) where T<:RestrictedTwoBodyOrbit = periapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

apoapsis_radius(orbit::T) where T<:RestrictedTwoBodyOrbit = apoapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

periapsis_velocity(orbit::T) where T<:RestrictedTwoBodyOrbit = speed(periapsis_radius(orbit), semimajor_axis(orbit), massparameter(system(orbit)))

apoapsis_velocity(orbit::T) where T<:RestrictedTwoBodyOrbit = speed(apoapsis_radius(orbit), semimajor_axis(orbit), massparameter(system(orbit)))

"""
Returns mass `m`.
"""
mass(body::RestrictedTwoBodySystem) = massparameter(body) / G

"""
Returns mass parameter `μ`.
"""
massparameter(body::RestrictedTwoBodySystem) = body.μ * massparameterunit(body)

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
period(orbit::RestrictedTwoBodyOrbit{C}) where C<:Union{Elliptical, Circular} = period(semimajor_axis(orbit), massparameter(system(orbit)))

true_anomoly(orbit::KeplerianOrbit) = orbit.state.ν * angularunit(orbit.state)

true_anomoly(orbit::CartesianOrbit) = true_anomoly(distance(orbit), specific_angular_momentum(orbit), eccentricity(orbit), massparameter(system(orbit)))

mean_motion(orbit::T) where T<:RestrictedTwoBodyOrbit = mean_motion(semimajor_axis(orbit), massparameter(system(orbit)))

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

specific_potential_energy(orbit::CartesianOrbit) = specific_potential_energy(distance(orbit.state), massparameter(system(orbit)))
specific_potential_energy(orbit::KeplerianOrbit) = specific_potential_energy(CartesianOrbit(orbit))

"""
Returns true if all elements are within `atol` of the other.
"""
function Base.isapprox(b1::RestrictedTwoBodySystem, b2::RestrictedTwoBodySystem; atol=1e-6)
    return ustrip(upreferred(abs(massparameter(b1) - massparameter(b2)))) < atol
end

"""
Returns true if all elements are identically equal.
"""
function Base.isequal(b1::RestrictedTwoBodySystem, b2::RestrictedTwoBodySystem)
    return massparameter(b1) == massparameter(b2)
end


=#