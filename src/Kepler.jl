#
# Restricted Two-body Problem models
#

"""
A mutable vector, with labels, for 6DOF Keplerian states.
"""
Base.@kwdef mutable struct OrbitalElements{F} <: AstrodynamicalState{F,6}
    e::F = 0.0
    a::F = 0.0
    i::F = 0.0
    Ω::F = 0.0
    ω::F = 0.0
    ν::F = 0.0

    OrbitalElements{F}(::UndefInitializer) where {F} = new{F}()
    OrbitalElements(::UndefInitializer) = OrbitalElements{Float64}(undef)

    OrbitalElements{F}(e, a, i, Ω, ω, ν) where {F} = new{F}(e, a, i, Ω, ω, ν)
    OrbitalElements(e, a, i, Ω, ω, ν) = new{promote_type(typeof(e), typeof(a), typeof(i), typeof(Ω), typeof(ω), typeof(ν))}(e, a, i, Ω, ω, ν)
    OrbitalElements{F}(state::NamedTuple) where {F} =
        let
            (; e, a, i, Ω, ω, ν) = merge((; e=zero(F), a=zero(F), i=zero(F), Ω=zero(F), ω=zero(F), ν=zero(F)), state)
            OrbitalElements{F}(e, a, i, Ω, ω, ν)
        end
    OrbitalElements(state::NamedTuple) = OrbitalElements{Float64}(state)
end

@doc OrbitalElements
const KeplerianState = OrbitalElements

"""
A parameter vector for Keplerian dynamics.
"""
Base.@kwdef struct KeplerianParameters{F} <: AstrodynamicalParameters{F,1}
    μ::F

    KeplerianParameters(μ) = new{typeof(μ)}(μ)
    KeplerianParameters{F}(μ) where {F} = new{F}(μ)
    KeplerianParameters(p::KeplerianParameters) = KeplerianParameters(p.μ)
    KeplerianParameters{F}(p::KeplerianParameters) where {F} = KeplerianParameters{F}(p.μ)
end

Base.convert(::Type{R2BParameters}, kepler::KeplerianParameters{B}) where {B} = R2BParameters{B}(kepler.μ)
Base.convert(::Type{R2BParameters{A}}, kepler::KeplerianParameters{B}) where {A,B} = R2BParameters{A}(kepler.μ)
Base.convert(::Type{KeplerianParameters}, kepler::R2BParameters{B}) where {B} = KeplerianParameters{B}(kepler.μ)
Base.convert(::Type{KeplerianParameters{A}}, kepler::R2BParameters{B}) where {A,B} = KeplerianParameters{A}(kepler.μ)

R2BParameters(k::KeplerianParameters) = R2BParameters{eltype(k)}(k.μ)
KeplerianParameters(r::R2BParameters) = KeplerianParameters{eltype(r)}(r.μ)

Base.@pure paradigm(::KeplerianParameters) = "Idealized Keplerian Dynamics"

"""
An `Orbit` which exists within Keplerian dynamics.
"""
const KeplerianOrbit = Orbit{<:KeplerianState,<:KeplerianParameters}
AstrodynamicalModels.KeplerianOrbit(state::AbstractVector, parameters::AbstractVector) = Orbit(state isa AstrodynamicalState ? state : KeplerianState(state), KeplerianParameters(parameters))
AstrodynamicalModels.KeplerianOrbit(; state::AbstractVector, parameters::AbstractVector) = Orbit(state isa AstrodynamicalState ? state : KeplerianState(state), KeplerianParameters(parameters))
