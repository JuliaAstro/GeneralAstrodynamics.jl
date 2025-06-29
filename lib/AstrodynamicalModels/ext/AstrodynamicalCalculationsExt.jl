module AstrodynamicalCalculationsExt

using AstrodynamicalModels
using AstrodynamicalCalculations
import AstrodynamicalCalculations as AC

gm(μ) = μ
gm(μ::Union{<:R2BParameters,<:KeplerianParameters}) = μ.μ

function AC.keplerian_to_cartesian(state::KeplerianState, μ)
    p = μ isa Number ? R2BParameters(μ) : μ
    u = AC.keplerian_to_cartesian(state.e, state.a, state.i, state.Ω, state.ω, state.ν, p.μ)

    return CartesianState(u)
end
AC.keplerian_to_cartesian(orbit::CartesianOrbit) =
    Orbit(AC.keplerian_to_cartesian(orbit.state, orbit.parameters), orbit.parameters)

function AC.cartesian_to_keplerian(state::CartesianState, μ)
    p = μ isa Number ? R2BParameters(μ) : μ
    u = AC.cartesian_to_keplerian(
        state.x,
        state.y,
        state.z,
        state.ẋ,
        state.ẏ,
        state.ż,
        p.μ,
    )

    return KeplerianState(u)
end
AC.cartesian_to_keplerian(orbit::CartesianOrbit) =
    Orbit(AC.cartesian_to_keplerian(orbit.state, orbit.parameters), orbit.parameters)

AC.eccentricity(state::KeplerianState) = state.e
AC.eccentricity(state::KeplerianState, μ) = AC.eccentricity(state)
AC.eccentricity(state::CartesianState, μ) =
    AC.eccentricity(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.eccentricity(orbit::Orbit) = AC.eccentricity(orbit.state, orbit.parameters)

AC.eccentricity_vector(state::CartesianState, μ) =
    AC.eccentricity_vector(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.eccentricity_vector(state::KeplerianState, μ) =
    AC.eccentricity_vector(keplerian_to_cartesian(state, μ), μ)
AC.eccentricity_vector(orbit::Orbit) = AC.eccentricity_vector(orbit.state, orbit.parameters)

AC.semimajor_axis(state::KeplerianState) = state.a
AC.semimajor_axis(state::KeplerianState, μ) = AC.semimajor_axis(state)
AC.semimajor_axis(state::CartesianState, μ) =
    AC.semimajor_axis(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.semimajor_axis(orbit::Orbit) = AC.semimajor_axis(orbit.state, orbit.parameters)

AC.inclination(state::KeplerianState) = state.i
AC.inclination(state::KeplerianState, μ) = AC.inclination(state)
AC.inclination(state::CartesianState, μ) = AC.inclination(cartesian_to_keplerian(state, μ))
AC.inclination(orbit::Orbit) = AC.inclination(orbit.state, orbit.parameters)

AC.right_ascension_ascending_node(state::KeplerianState) = state.Ω
AC.right_ascension_ascending_node(state::KeplerianState, μ) =
    AC.right_ascension_ascending_node(state)
AC.right_ascension_ascending_node(state::CartesianState, μ) =
    AC.right_ascension_ascending_node(cartesian_to_keplerian(state, μ), gm(μ))
AC.right_ascension_ascending_node(orbit::Orbit) =
    AC.right_ascension_ascending_node(orbit.state, orbit.parameters)

AC.argument_of_periapsis(state::KeplerianState) = state.ω
AC.argument_of_periapsis(state::KeplerianState, μ) = AC.argument_of_periapsis(state)
AC.argument_of_periapsis(state::CartesianState, μ) =
    AC.argument_of_periapsis(cartesian_to_keplerian(state, μ), gm(μ))
AC.argument_of_periapsis(orbit::Orbit) =
    AC.argument_of_periapsis(orbit.state, orbit.parameters)

AC.true_anomaly(state::KeplerianState) = state.ω
AC.true_anomaly(state::KeplerianState, μ) = AC.true_anomaly(state)
AC.true_anomaly(state::CartesianState, μ) =
    AC.true_anomaly(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.true_anomaly(orbit::Orbit) = AC.true_anomaly(orbit.state, orbit.parameters)

AC.semi_parameter(state::KeplerianState) = AC.semi_parameter(state.a, state.e)
AC.semi_parameter(state::KeplerianState, μ) = AC.semi_parameter(state)
AC.semi_parameter(state::CartesianState, μ) =
    AC.semi_parameter(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.semi_parameter(orbit::Orbit) = AC.semi_parameter(orbit.state, orbit.parameters)

AC.specific_energy(state::KeplerianState, μ) = AC.specific_energy(state.a, gm(μ))
AC.specific_energy(state::CartesianState, μ) =
    AC.specific_energy(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.specific_energy(orbit::Orbit) = AC.specific_energy(orbit.state, orbit.parameters)

AC.specific_potential_energy(state::KeplerianState) =
    AC.orbital_radius(semi_parameter(state.a, state.e), state.e, state.ν)
AC.specific_potential_energy(state::KeplerianState, μ) = AC.specific_potential_energy(state)
AC.specific_potential_energy(state::CartesianState, μ) = AC.specific_potential_energy(
    state.x,
    state.y,
    state.z,
    state.ẋ,
    state.ẏ,
    state.ż,
    gm(μ),
)
AC.specific_potential_energy(orbit::Orbit) =
    AC.specific_potential_energy(orbit.state, orbit.parameters)

AC.specific_angular_momentum(state::CartesianState, μ) = AC.specific_angular_momentum(
    state.x,
    state.y,
    state.z,
    state.ẋ,
    state.ẏ,
    state.ż,
    gm(μ),
)
AC.specific_angular_momentum(state::KeplerianState, μ) =
    AC.specific_angular_momentum(keplerian_to_cartesian(state, μ), μ)
AC.specific_angular_momentum(orbit::Orbit) =
    AC.specific_angular_momentum(orbit.state, orbit.parameters)

AC.specific_angular_momentum_vector(state::CartesianState, μ) =
    AC.specific_angular_momentum_vector(
        state.x,
        state.y,
        state.z,
        state.ẋ,
        state.ẏ,
        state.ż,
        gm(μ),
    )
AC.specific_angular_momentum_vector(state::KeplerianState, μ) =
    AC.specific_angular_momentum_vector(keplerian_to_cartesian(state, μ), μ)
AC.specific_angular_momentum_vector(orbit::Orbit) =
    AC.specific_angular_momentum_vector(orbit.state, orbit.parameters)

AC.c3(state::CartesianState, μ) =
    AC.c3(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.c3(state::KeplerianState, μ) = AC.c3(keplerian_to_cartesian(state, μ), μ)
AC.c3(orbit::Orbit) = AC.c3(orbit.state, orbit.parameters)

AC.v_infinity(state::CartesianState, μ) =
    AC.v_infinity(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.v_infinity(state::KeplerianState, μ) = AC.v_infinity(keplerian_to_cartesian(state, μ), μ)
AC.v_infinity(orbit::Orbit) = AC.v_infinity(orbit.state, orbit.parameters)

AC.orbital_period(state::KeplerianState, μ) = AC.orbital_period(state.a, gm(μ))
AC.orbital_period(state::CartesianState, μ) =
    AC.orbital_period(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.orbital_period(orbit::Orbit) = AC.orbital_period(orbit.state, orbit.parameters)

AC.orbital_radius(state::KeplerianState, μ) =
    AC.orbital_radius(AC.semi_parameter(state.a, state.e), state.e, state.ν)
AC.orbital_radius(state::CartesianState, μ) =
    AC.orbital_radius(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.orbital_radius(orbit::Orbit) = AC.orbital_radius(orbit.state, orbit.parameters)

AC.orbital_speed(state::KeplerianState, μ) =
    AC.orbital_speed(AC.orbital_radius(state, μ), state.a, gm(μ))
AC.orbital_speed(state::CartesianState, μ) =
    AC.orbital_speed(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.orbital_speed(orbit::Orbit) = AC.orbital_speed(orbit.state, orbit.parameters)

AC.mean_motion(state::KeplerianState, μ) = AC.mean_motion(state.a, gm(μ))
AC.mean_motion(state::CartesianState, μ) =
    AC.mean_motion(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.mean_motion(orbit::Orbit) = AC.mean_motion(orbit.state, orbit.parameters)

AC.apoapsis_radius(state::KeplerianState) = AC.apoapsis_radius(state.a, state.e)
AC.apoapsis_radius(state::KeplerianState, μ) = AC.apoapsis_radius(state)
AC.apoapsis_radius(state::CartesianState, μ) =
    AC.apoapsis_radius(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.apoapsis_radius(orbit::Orbit) = AC.apoapsis_radius(orbit.state, orbit.parameters)

AC.periapsis_radius(state::KeplerianState) = AC.periapsis_radius(state.a, state.e)
AC.periapsis_radius(state::KeplerianState, μ) = AC.periapsis_radius(state)
AC.periapsis_radius(state::CartesianState, μ) =
    AC.periapsis_radius(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, gm(μ))
AC.periapsis_radius(orbit::Orbit) = AC.periapsis_radius(orbit.state, orbit.parameters)

AC.conic(state::KeplerianState) = AC.conic(state.e)
AC.conic(state::KeplerianState, μ) = AC.conic(state)
AC.conic(state::CartesianState, μ) = AC.conic(AC.eccentricity(state, μ))
AC.conic(orbit::Orbit) = AC.conic(orbit.state, orbit.parameters)

function Base.show(io::IO, ::MIME"text/plain", orbit::Union{<:R2BOrbit,<:KeplerianOrbit})
    println(
        io,
        "$(AC.conic(orbit)) Orbit in $(AstrodynamicalModels.paradigm(orbit.parameters))\n",
    )

    for line in eachline(IOBuffer(string(orbit.state)))
        println(io, "  ", line)
    end

    println(io)

    for line in eachline(IOBuffer(string(orbit.parameters)))
        println(io, "  ", line)
    end

end

end
