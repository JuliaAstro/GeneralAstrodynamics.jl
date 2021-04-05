#
# Types for celestial bodies
#


NAIF_IDS = Dict{String, Int}(
    "Solar System Barycenter" => 0,
    "Sun" => 10,
    "Mercury Barycenter" => 1,
    "Mercury" => 199, 
    "Venus Barycenter" => 2,
    "Venus" => 299,
    "Earth-Moon Barycenter" => 3,
    "Moon" => 301,
    "Earth" => 399,
    "Mars Barycenter" => 4,
    "Mars" => 499,
     "Jupiter Barycenter" => 5,
     "Jupiter" => 599,
     "Saturn Barycenter" => 6,
     "Saturn" => 699,
     "Uranus Barycenter" => 7,
     "Uranus" => 799,
     "Neptune Barycenter" => 8,
     "Neptune" => 899,
     "Pluto Barycenter" => 9,
     "Pluto" => 999,
)

"""
Type representing large bodies in space. Currently, the following
solar system bodies are supported:

Sun, Mercury, Venus, Earth, Moon (Luna), Mars, Jupiter, 
Saturn, Uranus, Neptune, Pluto.
"""
struct CelestialBody{F<:AbstractFloat}
    r::SVector{3, Unitful.Length{F}}
    v::SVector{3, Unitful.Velocity{F}}
    mean_radius::Length{F}
    μ::MassParameter{F}
    name::String

    function CelestialBody(r::R, v::V, m::Mass{<:AbstractFloat}, radius::Length{<:AbstractFloat}, name::String="") where {
            R <: AbstractVecOrMat{<:Unitful.Length},
            V <: AbstractVecOrMat{<:Unitful.Velocity}
        }
        T = promote_type(typeof.(ustrip.(vcat(r,v)))..., typeof(ustrip(m)), typeof(ustrip(radius)))
        return new{T}(SVector{3,Unitful.Length{T}}(r), SVector{3,Unitful.Velocity{T}}(v), T(radius), T(G * m), name)
    end

    function CelestialBody(r::R, v::V, μ::MassParameter{<:AbstractFloat}, radius::Length{<:AbstractFloat}, name::String="") where {
            R <: AbstractVecOrMat{<:Unitful.Length},
            V <: AbstractVecOrMat{<:Unitful.Velocity}
        }
        T = promote_type(typeof.(ustrip.(vcat(r,v)))..., typeof(ustrip(μ)), typeof(ustrip(radius)))
        new{T}(SVector{3,Unitful.Length{T}}(r), SVector{3,Unitful.Velocity{T}}(v), radius, μ, name)
    end

    function CelestialBody(r::AbstractVecOrMat, v::AbstractVecOrMat, v::V, μ::T1, radius::T2, name::String="") where {T1<:AbstractFloat, T2<:AbstractFloat}
        @warn "No units provided! Assuming km and km^3/s^2."
        T = promote_type(T1, T2)
        return new{T}(T(radius), T(μ), name)
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

    CelestialBody(body::CelestialBody) = CelestialBody(body.μ, body.radius, body.name)

end

Base.convert(::Type{T}, b::CelestialBody) where {T<:AbstractFloat} = CelestialBody(T(b.μ), T(b.radius), b.name)
Base.promote(::Type{CelestialBody{A}}, ::Type{CelestialBody{B}}) where {A<:AbstractFloat, B<:AbstractFloat} = CelestialBody{promote_type(A,B)}
Core.Float16(o::CelestialBody) = convert(Float16, o)
Core.Float32(o::CelestialBody) = convert(Float32, o)
Core.Float64(o::CelestialBody) = convert(Float64, o)
Base.MPFradius.BigFloat(o::CelestialBody) = convert(BigFloat, o)

"""
Custom display for `CelestialBody` instances.
"""
function Base.show(io::IO, body::CelestialBody)

    println(io, crayon"blue", "CelestialBody:")
    println(io, crayon"default", 
                "    Mass:           ", ustrip(u"kg", body.μ / G), " ", u"kg")
    println(io, "    radiusadius:         ", ustrip(u"km", body.radius), " ", u"km")
    println(io, "    Mass Parameter: ", ustrip(u"km^3/s^2", body.μ), " ", u"km^3/s^2")

end
