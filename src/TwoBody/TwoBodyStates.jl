#
#   TwoBodyStates.jl
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
Abstract type for all two-body orbital representations.
"""
abstract type RestrictedTwoBodySystem{C<:AbstractConic, F<:AbstractFloat} <: OrbitalSystem end


"""
Type representing large bodies in space. Currently, the following
solar system bodies are supported:

Sun, Mercury, Venus, Earth, Moon (Luna), Mars, Jupiter, 
Saturn, Uranus, Neptune, Pluto.
"""
struct CelestialBody{F<:AbstractFloat}
    R::Length{F}
    μ::MassParameter{F}
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

Base.convert(::Type{T}, b::CelestialBody) where {T<:AbstractFloat} = CelestialBody(T(b.μ), T(b.R), b.name)
Base.promote(::Type{CelestialBody{A}}, ::Type{CelestialBody{B}}) where {A<:AbstractFloat, B<:AbstractFloat} = CelestialBody{promote_type(A,B)}
Core.Float16(o::CelestialBody) = convert(Float16, o)
Core.Float32(o::CelestialBody) = convert(Float32, o)
Core.Float64(o::CelestialBody) = convert(Float64, o)
Base.MPFR.BigFloat(o::CelestialBody) = convert(BigFloat, o)

"""
Custom display for `CelestialBody` instances.
"""
function Base.show(io::IO, body::CelestialBody)

    println(io, crayon"blue", "CelestialBody:")
    println(io, crayon"default", 
                "    Mass:           ", ustrip(u"kg", body.μ / G), " ", u"kg")
    println(io, "    Radius:         ", ustrip(u"km", body.R), " ", u"km")
    println(io, "    Mass Parameter: ", ustrip(u"km^3/s^2", body.μ), " ", u"km^3/s^2")

end

"""
Struct for storing `TwoBody` Cartesian states for all conics.
"""
struct TwoBodyState{C<:AbstractConic, F<:AbstractFloat} <: RestrictedTwoBodySystem{C,F}
    r::SVector{3, Length{F}}
    v::SVector{3, Velocity{F}}
    body::CelestialBody{F}

    function TwoBodyState(r::R, v::V, body::CelestialBody) where {R <: AbstractVector{<:Length}, V <: AbstractVector{<:Velocity}}
        C = conic(eccentricity(r, v, body.μ))
        T = promote_type(typeof(ustrip(r[1])), typeof(ustrip(v[1])), typeof(ustrip(body.μ)))
        if !(T <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            T = Float64
        end
        return new{C,T}(SVector{3,Length{T}}(r...), SVector{3, Velocity{T}}(v...), T(body))
    end

    function TwoBodyState(r::R, v::V, body::CelestialBody) where {R <: AbstractVector{<:Real}, V <: AbstractVector{<:Real}}
        C = conic(eccentricity(r, v, body.μ))
        T = promote_type(typeof(ustrip(r[1])), typeof(ustrip(v[1])), typeof(ustrip(body.μ)))
        if !(T <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            T = Float64
        end
        return new{C,T}(SVector{3,Length{T}}((r * u"km")...), SVector{3, Velocity{T}}((v * u"km/s")...), T(body))
    end

    TwoBodyState(r, v, μ::T) where T <: Real = TwoBodyState(r, v, CelestialBody(μ))
    TwoBodyState(orbit::TwoBodyState) = TwoBodyState(orbit.r, orbit.v, orbit.body)
end

Base.convert(::Type{T}, o::TwoBodyState) where {T<:AbstractFloat} = TwoBodyState(T.(o.r), T.(o.v), convert(T, o.body))
Base.promote(::Type{TwoBodyState{A}}, ::Type{TwoBodyState{B}}) where {A<:AbstractFloat, B<:AbstractFloat} = Orbit{promote_type(A,B)}
Core.Float16(o::TwoBodyState) = convert(Float16, o)
Core.Float32(o::TwoBodyState) = convert(Float32, o)
Core.Float64(o::TwoBodyState) = convert(Float64, o)
Base.MPFR.BigFloat(o::TwoBodyState) = convert(BigFloat, o)

"""
Alias for `TwoBodyState`.
"""
Orbit(r, v, body) = TwoBodyState(r, v, body)
Orbit(e, a, i, Ω, ω, ν, body) = KeplerianState(e, a, i, Ω, ω, ν, body)

"""
Struct for storing `Keplerian` states for all conics.
"""
struct KeplerianState{C<:AbstractConic, F<:AbstractFloat} <: RestrictedTwoBodySystem{C,F}
    e::F
    a::Length{F}
    i::DimensionlessQuantity{F}
    Ω::DimensionlessQuantity{F}
    ω::DimensionlessQuantity{F}
    ν::DimensionlessQuantity{F}
    body::CelestialBody{F}
    
    function KeplerianState(e::E, a::Length{A}, i::DimensionlessQuantity{I}, 
                            Ω::DimensionlessQuantity{O}, ω::DimensionlessQuantity{W}, 
                            ν::DimensionlessQuantity{V}, body::CelestialBody{B}) where {
            E <: Real,
            A <: Real,
            I <: Real,
            O <: Real,
            W <: Real,
            V <: Real,
            B <: Real
        }
        C = conic(e)
        T = promote_type(typeof(e), typeof(ustrip(a)), typeof(ustrip(i)), typeof(ustrip(Ω)), 
                         typeof(ustrip(ω)), typeof(ustrip(ν)), typeof(ustrip(body.μ)))
        if !(T <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            T = Float64
        end
        return new{C,T}(T(e), T(a), T(i), T(Ω), T(ω), T(ν), T(body))
    end

    function KeplerianState(e::E, a::A, i::I, 
                            Ω::O, ω::W, ν::V,
                            body::CelestialBody{B}) where {
            E <: Real,
            A <: Real,
            I <: Real,
            O <: Real,
            W <: Real,
            V <: Real,
            B <: Real
        }
        @warn "No units provided! Assuming km and rad."
        return KeplerianState(e, a * u"km", (i % 2π) * u"rad", (Ω % 2π) * u"rad", 
                              (ω % 2π) * u"rad", (ν % 2π) * u"rad", body)
    end

    KeplerianState(e, a, i, Ω, ω, ν, μ::T) where T <: Real = KeplerianState(e, a, i, Ω, ω, ν, CelestialBody(μ))
    KeplerianState(orbit::KeplerianState) = KeplerianState(orbit.e, orbit.a, orbit.i, orbit.Ω, orbit.ω, orbit.ν, orbit.body)
end

Base.convert(::Type{T}, o::KeplerianState) where {T<:AbstractFloat} = KeplerianState(T(o.e), T(o.a), T(o.i), T(o.Ω), T(o.ω), T(o.ν), convert(T, o.body))
Base.promote(::Type{KeplerianState{A}}, ::Type{KeplerianState{B}}) where {A<:AbstractFloat, B<:AbstractFloat} = KeplerianState{promote_type(A,B)}
Core.Float16(o::KeplerianState) = convert(Float16, o)
Core.Float32(o::KeplerianState) = convert(Float32, o)
Core.Float64(o::KeplerianState) = convert(Float64, o)
Base.MPFR.BigFloat(o::KeplerianState) = convert(BigFloat, o)

TwoBodyState(orbit::KeplerianState) = TwoBodyState(cartesian(orbit)..., orbit.body)
KeplerianState(orbit::TwoBodyState) = KeplerianState(keplerian(orbit)..., orbit.body)

"""
Custom display for TwoBodyState instances.
"""
function Base.show(io::IO, orbit::TwoBodyState)

    println(io, conic(orbit), " Two-body State (Cartesian):")
    println(io, "")

    println(io, "    Position (inertial):    [", 
                ustrip(u"km", orbit.r[1]), ", ", 
                ustrip(u"km", orbit.r[2]), ", ", 
                ustrip(u"km", orbit.r[3]), "] ", u"km")
    println(io, "    Velocity (inertial):    [", 
                ustrip(u"km/s", orbit.v[1]), ", ", 
                ustrip(u"km/s", orbit.v[2]), ", ", 
                ustrip(u"km/s", orbit.v[3]), "] ", u"km/s")

    println(io, "")
    println(io, "    $(orbit.body.name == "" ? "Body" : orbit.body.name) (μ):               ",
                ustrip(u"km^3 / s^2", orbit.body.μ), " ", u"km^3/s^2")
end

"""
Custom display for KeplerianState instances.
"""
function Base.show(io::IO, orbit::KeplerianState)

    println(io, conic(orbit), " Two-body State (Keplerian):")
    println(io, "")
    
    println(io, "")
    println(io, "    Eccentricity:           ", 
                orbit.e)    
    println(io, "    Semimajor Axis:         ", 
                ustrip(u"km", orbit.a), " ", u"km")
    println(io, "    Inclination:            ", 
                ustrip(u"°", orbit.i), u"°")
    println(io, "    RAAN:                   ", 
                ustrip(u"°", orbit.Ω), u"°")
    println(io, "    Arg. Periapsis:         ", 
                ustrip(u"°", orbit.ω), u"°")
    println(io, "    True Anomoly:           ", 
                ustrip(u"°", orbit.ν), u"°")

    println(io, "")
    println(io, "    $(orbit.body.name == "" ? "Body" : orbit.body.name) (μ):               ",
                ustrip(u"km^3 / s^2", orbit.body.μ), " ", u"km^3/s^2")
end

# Constants

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
    
