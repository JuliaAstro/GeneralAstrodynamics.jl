#
# Methods for orbit states and state vectors.
#

Calculations.keplerian(state::CartesianStateVector, μ) = keplerian(position(state), velocity(state), μ)
Calculations.keplerian(orbit::CartesianR2BPOrbit) = keplerian(state(orbit), massparameter(system(orbit)))
Calculations.keplerian(orbit::KeplerianR2BPOrbit) = orbit

Calculations.cartesian(orbit::CartesianR2BPOrbit) = orbit
Calculations.cartesian(orbit::KeplerianR2BPOrbit) = cartesian(state(orbit), massparameter(system(orbit)))
Calculations.cartesian(state::KeplerianState, μ) = cartesian(
    eccentricity(state),
    semimajor_axis(state),
    inclination(state),
    RAAN(state),
    argument_of_periapsis(state),
    true_anomoly(state),
    μ
)

Base.position(state::CartesianStateVector) = States.get_r(state) * lengthunit(state)
Base.position(orbit::KeplerianR2BPOrbit) = position(CartesianState(cartesian(state(orbit), massparameter(state(orbit)))))
Base.position(orbit::Orbit) = position(state(orbit))

Calculations.distance(state::CartesianStateVector) = norm(position(state))
Calculations.distance(orbit::Orbit) = norm(position(orbit))

velocity(state::CartesianStateVector) = States.get_v(state) * velocityunit(state)
velocity(orbit::Orbit) = velocity(state(orbit))

Calculations.speed(state::CartesianStateVector) = norm(velocity(state))
Calculations.speed(orbit::KeplerianR2BPOrbit) = speed(semi_parameter(semimajor_axis(state(orbit)), eccentricity(state(orbit))), semimajor_axis(state(orbit)), massparameter(system(orbit)))
Calculations.speed(orbit::Orbit) = norm(velocity(state(orbit)))

Calculations.eccentricity(state::KeplerianState) = States.get_e(state)
Calculations.eccentricity(orbit::KeplerianR2BPOrbit) = eccentricity(state(orbit))
Calculations.eccentricity(orbit::CartesianR2BPOrbit) = eccentricity(position(state(orbit)), velocity(state(orbit)), massparameter(system(orbit)))

Calculations.semimajor_axis(state::KeplerianState) = States.get_a(state) * lengthunit(state)
Calculations.semimajor_axis(orbit::CartesianR2BPOrbit) = semimajor_axis(distance(state(orbit)), speed(state(orbit)), massparameter(system(orbit)))
Calculations.semimajor_axis(orbit::KeplerianR2BPOrbit) = semimajor_axis(state(orbit))

inclination(state::KeplerianState) = States.get_i(state) * angularunit(state)
inclination(orbit::KeplerianR2BPOrbit) = inclination(state(orbit))

RAAN(state::KeplerianState) = States.get_Ω(state) * angularunit(state)
RAAN(orbit::KeplerianR2BPOrbit) = RAAN(state(orbit))

argument_of_periapsis(state::KeplerianState) = States.get_ω(state) * angularunit(state)
argument_of_periapsis(orbit::KeplerianR2BPOrbit) = argument_of_periapsis(state(orbit))

Calculations.true_anomoly(state::KeplerianState) = States.get_ν(state) * angularunit(state)
Calculations.true_anomoly(orbit::KeplerianR2BPOrbit) = true_anomoly(state(orbit))
Calculations.true_anomoly(orbit::CartesianR2BPOrbit) = true_anomoly(distance(orbit), specific_angular_momentum(orbit), eccentricity(orbit), massparameter(system(orbit)))

Calculations.conic(orbit::R2BPOrbit) = conic(eccentricity(orbit))

Calculations.specific_angular_momentum_vector(orbit::CartesianR2BPOrbit) = specific_angular_momentum_vector(position(state(orbit)), velocity(state(orbit)))

Calculations.specific_angular_momentum(orbit::CartesianR2BPOrbit) = specific_angular_momentum(position(state(orbit)), velocity(state(orbit)))

Calculations.specific_energy(orbit::CartesianR2BPOrbit) = specific_energy(distance(state(orbit)), speed(state(orbit)), massparameter(system(orbit)))
Calculations.specific_energy(orbit::KeplerianR2BPOrbit) = specific_energy(semimajor_axis(orbit), massparameter(system(orbit)))

Calculations.C3(orbit::R2BPOrbit) = C3(distance(orbit), speed(orbit), massparameter(system(orbit)))

Calculations.v_infinity(orbit::R2BPOrbit) = v_infinity(distance(orbit), speed(orbit), massparameter(system(orbit)))

Calculations.eccentricity_vector(orbit::CartesianR2BPOrbit) = eccentricity_vector(position(state(orbit)), velocity(state(orbit)), massparameter(system(orbit)))

Calculations.semi_parameter(orbit::R2BPOrbit) = semi_parameter(semimajor_axis(orbit), eccentricity(orbit))

Calculations.periapsis_radius(orbit::R2BPOrbit) = periapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

Calculations.apoapsis_radius(orbit::R2BPOrbit) = apoapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

periapsis_velocity(orbit::R2BPOrbit) = speed(periapsis_radius(orbit), semimajor_axis(orbit), massparameter(system(orbit)))

apoapsis_velocity(orbit::R2BPOrbit) = speed(apoapsis_radius(orbit), semimajor_axis(orbit), massparameter(system(orbit)))

Calculations.period(orbit::R2BPOrbit) = conic(orbit) ∈ (Circular, Elliptical) ? period(semimajor_axis(orbit), massparameter(system(orbit))) : conic(orbit) == Parabolic ? Inf * timeunit(orbit) : NaN * timeunit(orbit)

Calculations.mean_motion(orbit::R2BPOrbit) = mean_motion(semimajor_axis(orbit), massparameter(system(orbit)))

function mean_motion_vector(orbit::R2BPOrbit) 
#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3}(0, 0, 1)
    return k̂ × specific_angular_momentum_vector(orbit)
end

Calculations.time_since_periapsis(orbit::R2BPOrbit) = time_since_periapsis(mean_motion(orbit), eccentricity(orbit), eccentric_anomoly(orbit))

Calculations.specific_potential_energy(orbit::CartesianR2BPOrbit) = specific_potential_energy(distance(orbit.state), massparameter(system(orbit)))

Calculations.jacobi_constant(orbit::CR3BPOrbit) = jacobi_constant(States.get_r(state(orbit)), States.get_v(state(orbit)), massparameter(system(orbit)))

Calculations.potential_energy(orbit::CR3BPOrbit) = potential_energy(States.get_r(state(orbit)), massparameter(system(orbit)))

distance_to_primary(orbit::CR3BPOrbit) = norm(States.get_r(state(orbit)) .- SVector(-massparameter(orbit), 0, 0)) 

distance_to_secondary(orbit::CR3BPOrbit) = norm(States.get_r(state(orbit)) .- SVector(1-massparameter(orbit), 0, 0)) 

function Calculations.kepler(orbit::Orbit, Δt = period(orbit); kwargs...) 
    r, v = kepler(position(orbit), velocity(orbit), massparameter(orbit), Δt; kwargs...)
    cart = CartesianState(r, v)
    return Orbit(cart, system(orbit))
end