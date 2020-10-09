#
#   TwoBodyStates.jl
#
#   Describes Two Body Orbits through Cartesian coordinates and Orbital Elements.
# 

"""
    AbstractConic

Abstract type for all four conic sections.
"""
abstract type AbstractConic end


"""
    TwoBodyState

Abstract type for all two-body orbital representations.
"""
abstract type TwoBodySystem{T<:AbstractConic} <: OrbitalSystem end


"""
    Circular

Type for orbits in the circular conic section.
"""
struct Circular <: AbstractConic end

"""
    Elliptical

Type for orbits in the elliptical conic section.
"""
struct Elliptical <: AbstractConic end

"""
    Parabolic

Type for orbits in the parabolic conic section.
"""
struct Parabolic <: AbstractConic end

"""
    Hyperbolic

Type for orbits in the hyperbolic conic section.
"""
struct Hyperbolic <: AbstractConic end

"""
    Invalid

Type for invalid orbits (orbits with NaN fields)
"""
struct Invalid <: AbstractConic end

"""
    struct CelestialBody(m, R, μ)
    CelestialBody(m, R) = CelestialBody(m, R, G * m)

Type representing large bodies in space. Currently, the following
solar system bodies are supported:

Sun, Mercury, Venus, Earth, Moon (Luna), Mars, Jupiter, 
Saturn, Uranus, Neptune, Pluto.
"""
struct CelestialBody
    m::Quantity
    R::Quantity
    μ::Quantity
end
CelestialBody(m, R) = CelestialBody(m, R, G * m)


# All data pulled from the following references:
# [1] https://en.wikipedia.org/wiki/List_of_Solar_System_objects_by_size
# [2] https://docs.astropy.org/en/stable/constants/#module-astropy.constants
Sun = CelestialBody(
    1.98840987e30u"kg",
    696342.0u"km")
Mercury = CelestialBody(
    330.1e21u"kg",
    2439.7u"km")
Venus = CelestialBody(
    4867.5e21u"kg",
    6051.8u"km")
Earth = CelestialBody(
    5.97216787e24u"kg",
    6371.0u"km")
Moon = CelestialBody(
    73.42e21u"kg",
    1737.4u"km")
Luna = Moon
Mars = CelestialBody(
    641.7e21u"kg",
    3389.5u"km")
Jupiter = CelestialBody(
    1.8981246e27u"kg",
    69911.0u"km")
Saturn = CelestialBody(
    568340e21u"kg",
    58232.0u"km")
Uranus = CelestialBody(
    86813e21u"kg",
    25362.0u"km")
Neptune = CelestialBody(
    102413e21u"kg",
    24622.0u"km")
Pluto = CelestialBody(
    13.03e21u"kg",
    1188.3u"km")
    
"""
    Orbit{T<:AbstractConic}

Struct for storing TwoBody orbital states for all conics.
"""
struct Orbit{T<:AbstractConic} <: TwoBodySystem{T}

    # Cartesian representation
    r̅::SVector{3, Unitful.Length{Float64}}
    v̅::SVector{3, Unitful.Velocity{Float64}}

    # Keplerian representation
    e::Union{Float64,Unitful.DimensionlessQuantity{Float64}}
    a::Unitful.Length{Float64}
    i::Unitful.DimensionlessQuantity{Float64}
    Ω::Unitful.DimensionlessQuantity{Float64}
    ω::Unitful.DimensionlessQuantity{Float64}
    ν::Unitful.DimensionlessQuantity{Float64}

    # Body
    body::CelestialBody

end

"""
    InvalidOrbit(body::CelestialBody)

Returns a `Orbit` with `NaN` state values. Used by 
`propagate_twobody` and `kepler` to indicate failed convergance.
"""
InvalidOrbit(body::CelestialBody) = Orbit(SVector{3}(NaN * u"km", NaN * u"km", NaN * u"km"),
                                                        SVector{3}(NaN * u"km", NaN * u"km", NaN * u"km"),
                                                        body)

"""
    isinvalid(orbit::Orbit)

Checks for `NaN` valued orbital states, which are used to
indicate an invalid `Orbit`.
"""
isinvalid(orbit::Orbit)    = all(map(x->!isnan(getfield(orbit, x)), [:r̅, :v̅, :e, :a, :i, :Ω, :ω, :ν])) 