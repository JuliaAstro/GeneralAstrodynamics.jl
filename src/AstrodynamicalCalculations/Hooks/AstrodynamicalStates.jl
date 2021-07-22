#
# Defines AstrodynamicalCalculations methods for 
# AstrodynamicalStates types
#

keplerian(state::CartesianState, μ) = keplerian(position(state), velocity(state), μ)
keplerian(orbit::CartesianR2BPOrbit) = keplerian(state(orbit), massparameter(system(orbit)))
keplerian(orbit::KeplerianR2BPOrbit) = orbit

cartesian(orbit::CartesianR2BPOrbit) = orbit
cartesian(orbit::KeplerianR2BPOrbit) = cartesian(state(orbit), massparameter(system(orbit)))
cartesian(state::KeplerianState, μ) = cartesian(
    eccentricity(state),
    semimajor_axis(state),
    inclination(state),
    RAAN(state),
    argument_of_periapsis(state),
    true_anomoly(state),
    μ
)

Base.position(state::CartesianState) = AstrodynamicalStates.get_r(state) * lengthunit(state)
Base.position(orbit::CartesianR2BPOrbit) = position(state(orbit))
Base.position(orbit::KeplerianR2BPOrbit) = position(cartesian(state(orbit), massparameter(state(orbit))))

distance(state::CartesianState) = norm(position(state))
distance(orbit::Orbit) = norm(position(orbit))

velocity(state::CartesianState) = AstrodynamicalStates.get_v(state) * velocityunit(state)
velocity(orbit::CartesianR2BPOrbit) = velocity(state(orbit))

speed(state::CartesianState) = norm(velocity(state))
speed(orbit::CartesianR2BPOrbit) = norm(velocity(state(orbit)))
speed(orbit::KeplerianR2BPOrbit) = speed(semi_parameter(semimajor_axis(state(orbit)), eccentricity(state(orbit))), semimajor_axis(state(orbit)), massparameter(system(orbit)))

eccentricity(state::KeplerianState) = AstrodynamicalStates.get_e(state)
eccentricity(orbit::KeplerianR2BPOrbit) = eccentricity(state(orbit))
eccentricity(orbit::CartesianR2BPOrbit) = eccentricity(position(state(orbit)), velocity(state(orbit)), massparameter(system(orbit)))

semimajor_axis(state::KeplerianState) = AstrodynamicalStates.get_a(state) * lengthunit(state)
semimajor_axis(orbit::CartesianR2BPOrbit) = semimajor_axis(distance(state(orbit)), speed(state(orbit)), massparameter(system(orbit)))
semimajor_axis(orbit::KeplerianR2BPOrbit) = semimajor_axis(state(orbit))

inclination(state::KeplerianState) = AstrodynamicalStates.get_i(state) * angularunit(state)
inclination(orbit::KeplerianR2BPOrbit) = inclination(state(orbit))

RAAN(state::KeplerianState) = AstrodynamicalStates.get_Ω(state) * angularunit(state)
RAAN(orbit::KeplerianR2BPOrbit) = RAAN(state(orbit))

argument_of_periapsis(state::KeplerianState) = AstrodynamicalStates.get_ω(state) * angularunit(state)
argument_of_periapsis(orbit::KeplerianR2BPOrbit) = argument_of_periapsis(state(orbit))

true_anomoly(state::KeplerianState) = AstrodynamicalStates.get_ν(state) * angularunit(state)
true_anomoly(orbit::KeplerianR2BPOrbit) = true_anomoly(state(orbit))
true_anomoly(orbit::CartesianR2BPOrbit) = true_anomoly(distance(orbit), specific_angular_momentum(orbit), eccentricity(orbit), massparameter(system(orbit)))

conic(orbit::R2BPOrbit) = conic(eccentricity(orbit))

specific_angular_momentum_vector(orbit::CartesianR2BPOrbit) = specific_angular_momentum_vector(position(state(orbit)), velocity(state(orbit)))

specific_angular_momentum(orbit::CartesianR2BPOrbit) = specific_angular_momentum(position(state(orbit)), velocity(state(orbit)))

specific_energy(orbit::CartesianR2BPOrbit) = specific_energy(distance(state(orbit)), speed(state(orbit)), massparameter(system(orbit)))
specific_energy(orbit::KeplerianR2BPOrbit) = specific_energy(semimajor_axis(orbit), massparameter(system(orbit)))

C3(orbit::R2BPOrbit) = C3(distance(orbit), speed(orbit), massparameter(system(orbit)))

v_infinity(orbit::R2BPOrbit) = v_infinity(distance(orbit), speed(orbit), massparameter(system(orbit)))

eccentricity_vector(orbit::CartesianR2BPOrbit) = eccentricity_vector(position(state(orbit)), velocity(state(orbit)), massparameter(system(orbit)))

semi_parameter(orbit::R2BPOrbit) = semi_parameter(semimajor_axis(orbit), eccentricity(orbit))


periapsis_radius(orbit::R2BPOrbit) = periapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

apoapsis_radius(orbit::R2BPOrbit) = apoapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

periapsis_velocity(orbit::R2BPOrbit) = speed(periapsis_radius(orbit), semimajor_axis(orbit), massparameter(system(orbit)))

apoapsis_velocity(orbit::R2BPOrbit) = speed(apoapsis_radius(orbit), semimajor_axis(orbit), massparameter(system(orbit)))

period(orbit::R2BPOrbit) = conic(orbit) ∈ (Circular, Elliptical) ? period(semimajor_axis(orbit), massparameter(system(orbit))) : conic(orbit) == Parabolic ? Inf * timeunit(orbit) : NaN * timeunit(orbit)

mean_motion(orbit::R2BPOrbit) = mean_motion(semimajor_axis(orbit), massparameter(system(orbit)))

function mean_motion_vector(orbit::R2BPOrbit) 
#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3}(0, 0, 1)
    return k̂ × specific_angular_momentum_vector(orbit)
end

time_since_periapsis(orbit::R2BPOrbit) = time_since_periapsis(mean_motion(orbit), eccentricity(orbit), eccentric_anomoly(orbit))

specific_potential_energy(orbit::CartesianR2BPOrbit) = specific_potential_energy(distance(orbit.state), massparameter(system(orbit)))
