module AstrodynamicalCalculationsExt

using AstrodynamicalModels, AstrodynamicalCalculations

function AstrodynamicalCalculations.keplerian_to_cartesian(state::KeplerianState, μ)
    p = μ isa Number ? R2BParameters(μ) : μ
    u = keplerian_to_cartesian(state.e, state.a, state.i, state.Ω, state.ω, state.ν, p.μ)

    return CartesianState(u)
end

function AstrodynamicalCalculations.cartesian_to_keplerian(state::CartesianState, μ)
    p = μ isa Number ? R2BParameters(μ) : μ
    u = cartesian_to_keplerian(state.x, state.y, state.z, state.ẋ, state.ẏ, state.ż, p.μ)

    return KeplerianState(u)
end

Base.convert(::Type{CartesianOrbit}, orbit::KeplerianOrbit) = Orbit(keplerian_to_cartesian(orbit.state, orbit.parameters), orbit.parameters)
Base.convert(::Type{KeplerianOrbit}, orbit::CartesianOrbit) = Orbit(cartesian_to_keplerian(orbit.state, orbit.parameters), orbit.parameters)

AstrodynamicalModels.CartesianOrbit(orbit::KeplerianOrbit) = convert(CartesianOrbit, orbit)
AstrodynamicalModels.KeplerianOrbit(orbit::CartesianOrbit) = convert(KeplerianOrbit, orbit)
end