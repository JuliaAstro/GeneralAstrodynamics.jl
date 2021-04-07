#
#   RestrictedTwoBodySystems.jl
#
#   Describes Two Body Orbits through Cartesian coordinates and Orbital Elements.
# 

"""
Abstract type for all four conic sections.
"""
abstract type AbstractConic end

"""
Type for orbits in the circular conic section.
"""
struct Circular <: AbstractConic end

"""
Type for orbits in the elliptical conic section.
"""
struct Elliptical <: AbstractConic end

"""
Type for orbits in the parabolic conic section.
"""
struct Parabolic <: AbstractConic end

"""
Type for orbits in the hyperbolic conic section.
"""
struct Hyperbolic <: AbstractConic end

"""
Type for invalid orbits (orbits with NaN fields)
"""
struct Invalid <: AbstractConic end

"""
Type representing large bodies in space. Currently, the following
solar system bodies are supported:

Sun, Mercury, Venus, Earth, Moon (Luna), Mars, Jupiter, 
Saturn, Uranus, Neptune, Pluto.
"""
struct CelestialBody{F<:AbstractFloat, LU, MU} <: AbstractBody where {LU <: Unitful.LengthUnits, MU <: MassParameterUnits}
    R::F
    μ::F
    name::String

    function CelestialBody(m::Mass{<:AbstractFloat}, R::Length{<:AbstractFloat}, name::String="")
        T = promote_type(typeof(ustrip(m)), typeof(ustrip(R)))
        return new{T}(T(R), T(G * m), name)
    end

    function CelestialBody(μ::MassParameter{<:AbstractFloat}, R::Length{<:AbstractFloat}, name::String="")
        T = promote_type(typeof(ustrip(μ)), typeof(ustrip(R)))
        new{T}(R, μ, name)
    end

    function CelestialBody(μ::T1, R::T2, name::String="") where {T1<:AbstractFloat, T2<:AbstractFloat}
        @warn "No units provided! Assuming km and km^3/s^2."
        T = promote_type(T1, T2)
        return new{T}(T(R), T(μ), name)
    end


    function CelestialBody(μ::MassParameter, name::String="")
        @warn "No radius provided! Setting to NaN."
        return CelestialBody(μ, NaN * u"km", name)
    end

    function CelestialBody(μ::T, name::String="") where T<:AbstractFloat
        @warn "No units provided! Assuming km^3/s^2."
        return CelestialBody(μ * u"km^3/s^2", name)
    end

    CelestialBody(m::Mass) = CelestialBody(m * G)

    CelestialBody(body::CelestialBody) = CelestialBody(body.μ, body.R, body.name)

end

"""
Custom display for `CelestialBody` instances.
"""
function Base.show(io::IO, body::CelestialBody)

    println(io, "CelestialBody:")
    println(io, "    Mass:           ", ustrip(u"kg", body.μ / G), " ", u"kg")
    println(io, "    Radius:         ", ustrip(u"km", body.R), " ", u"km")
    println(io, "    Mass Parameter: ", ustrip(u"km^3/s^2", body.μ), " ", u"km^3/s^2")

end

"""
Abstract type for Keplerian states.
"""
abstract type AbstractKeplerianState{F<:AbstractFloat} <: FieldVector{6,F} end

"""
Struct for storing `Keplerian` states.
"""
struct KeplerianState{F<:AbstractFloat, LU, AU} <: AbstractKeplerianState{F} where {U, LU <: Unitful.LengthUnits, AU <: Unitful.Units{U, NoDims, nothing}}
    e::F
    a::F
    i::F
    Ω::F
    ω::F
    ν::F


    function KeplerianState(e::E, a::A, i::I, Ω::O, ω::W, ν::V; lengthunit = u"km", angularunit = u"rad") where {
        E <: Real, A <: Real, I <: Real, O <: Real, W <: Real, V <: Real
    }
        F = promote_type(E, A, I, O, W, V)
        if !(F <: AbstractFloat)
            @warn "Promoted type ($(string(F)) is not a float: defaulting to Float64."
            F = Float64
        end
        return new{F, lengthunit, angularunit}(e, a, i, Ω, ω, ν)
    end

    function KeplerianState(e::Real, a::Length, i::DimensionlessQuantity, 
                            Ω::DimensionlessQuantity, ω::DimensionlessQuantity, 
                            ν::DimensionlessQuantity)

        F = promote_type(typeof(e), ustrip(a), ustrip(i), ustrip(Ω), ustrip(ω), ustrip(ν))
        if !(F <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            F = Float64
        end
        lengthunit  = unit(a)
        angularunit = u"rad" ∈ unit.(i, Ω, ω, ν) ? u"rad" : unit(i)
        return KeplerianState(e, ustrip(lengthunit, a), ustrip(angularunit(i)), 
                              ustrip(angularunit(Ω)), ustrip(angularunit(ω)), ustrip(angularunit(ν));
                              lengthunit = lengthunit, angularunit = angularunit) 
    end

end

"""
Returns the `Unitful.Length` unit associated with the Keplerian state.
"""
Base.@pure lengthunit(::K) where K <: KeplerianState = C.parameters[2]

"""
Returns the dimmensionless unit associated with the Keplerian state.
"""
Base.@pure angularunit(::K) where K <: KeplerianState = C.parameters[3]

"""
An orbital state within the Restricted Two-body Problem.
"""
struct RestrictedTwoBodySystem{C<:AbstractConic, F<:AbstractFloat, T<:Union{CartesianState{F}, KeplerianState{F}}} <: AbstractOrbitalSystem
    state::T
    body::CelestialBody{F}

    function RestrictedTwoBodySystem(r, v, μ) 
        F = promote_type(eltype(r), eltype(v), typeof(μ))
        if !(F <: AbstractFloat)
            @warn "Promoted type $(string(F)) is not of type float. Defaulting to Float64."
            F = Float64
        end
        return RestrictedTwoBodySystem{conic(eccentricity(r,v,μ)), F, CartesianState{F}}(CartesianState(F.(r), F.(v)), CelestialBody(F(μ)))
    end

    function RestrictedTwoBodySystem(e, a, i, Ω, ω, ν, μ) 
        F = promote_type(typeof(e), typeof(a), typeof(i), typeof(Ω), typeof(ω), typeof(ν), typeof(μ))
        if !(F <: AbstractFloat)
            @warn "Promoted type $(string(F)) is not of type float. Defaulting to Float64."
            F = Float64
        end
        return RestrictedTwoBodySystem{conic(eccentricity(e)), F, KeplerianState{F}}(KeplerianState(F(e), F(a), F(i), F(Ω), F(ω), F(ν)), CelestialBody(F(μ)))
    end

end

"""
Alias for `RestrictedTwoBodySystem`.
"""
Orbit(r, v, body) = RestrictedTwoBodySystem(r, v, body)
Orbit(e, a, i, Ω, ω, ν, body) = KeplerianState(e, a, i, Ω, ω, ν, body)

"""
Custom display for KeplerianState instances.
"""
function Base.show(io::IO, orbit::KeplerianState)

    L = lengthunit(orbit)
    A = angularunit(orbit)

    println(io, conic(orbit), " Two-body State (Keplerian):")
    println(io, "")
    
    println(io, "")
    println(io, "    Eccentricity:           ", 
                orbit.e)    
    println(io, "    Semimajor Axis:         ", 
                orbit.a, " ", string(L))
    println(io, "    Inclination:            ", 
                orbit.i, A == u"rad" ? " " : "", string(A))
    println(io, "    RAAN:                   ", 
                orbit.Ω, A == u"rad" ? " " : "", string(A))
    println(io, "    Arg. Periapsis:         ", 
                orbit.ω, A == u"rad" ? " " : "", string(A))
    println(io, "    True Anomoly:           ", 
                orbit.ν, A == u"rad" ? " " : "", string(A))

end

# Constants

const KeplerianOrbit = RestrictedTwoBodySystem{C, F, KeplerianState{F}} where {C,F}
const CartesianOrbit = RestrictedTwoBodySystem{C, F, CartesianState{F}} where {C,F}

# All data pulled from the following references:
# [1] https://en.wikipedia.org/wiki/List_of_Solar_System_objects_by_size
# [2] https://docs.astropy.org/en/stable/constants/#module-astropy.constants

"""
Constant `CelestialBody` for our sun!
"""
const Sun = CelestialBody(1.327124400419393e11u"km^3/s^2", 696000.0u"km", "Sun")

"""
Constant `CelestialBody` for Mercury.
"""
const Mercury = CelestialBody(22031.78000000002u"km^3/s^2", 2439.7u"km", "Mercury")

"""
Constant `CelestialBody` for Venus.
"""
const Venus = CelestialBody(324858.592u"km^3/s^2", 6051.8u"km", "Venus")

"""
Constant `CelestialBody` for your home planet!
"""
const Earth = CelestialBody(398600.4354360959u"km^3/s^2", 6371.008366666666u"km", "Earth")

"""
Constant `CelestialBody` for our moon.
"""
const Moon = CelestialBody(4902.800066163796u"km^3/s^2", 1737.4000000000003u"km", "Moon")

"""
Constant `CelestialBody` (alias for our mooon).
"""
const Luna = Moon

"""
Constant `CelestialBody` for Mars.
"""
const Mars = CelestialBody(42828.37362069909u"km^3/s^2", 3389.5266666666666u"km", "Mars")

"""
Constant `CelestialBody` for Jupiter.
"""
const Jupiter = CelestialBody(1.2668653492180079e8u"km^3/s^2", 69946.0u"km", "Jupiter")

"""
Constant `CelestialBody` for Saturn.
"""
const Saturn = CelestialBody(3.793120749865224e7u"km^3/s^2", 58300.0u"km", "Saturn")

"""
Constant `CelestialBody` for Uranus.
"""
const Uranus = CelestialBody(5.793951322279009e6u"km^3/s^2", 25363.666666666668u"km", "Uranus")

"""
Constant `CelestialBody` for Neptune.
"""
const Neptune = CelestialBody(6.835099502439672e6u"km^3/s^2", 24623.0u"km", "Neptune")

"""
Constant `CelestialBody` for Pluto. We couldn't leave you out again!
"""
const Pluto = CelestialBody(869.6138177608749u"km^3/s^2", 1195.0u"km", "Pluto")
    
