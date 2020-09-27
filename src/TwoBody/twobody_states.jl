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
    CelestialBody(μ, R)

Type representing large bodies in space. Currently, only Earth 
and the Sun are supported. All bodies are treated as point 
masses. 

"""
struct CelestialBody
    μ::Quantity
    R::Quantity
end

Earth = CelestialBody(upreferred(1.0u"GMearth"), upreferred(1.0u"Rearth"))
Sun   = CelestialBody(upreferred(1.0u"GMsun"), upreferred(1.0u"Rsun"))

"""
    TwoBodyOrbit{T<:AbstractConic}

Struct for storing TwoBody orbital states for all conics.
"""
struct TwoBodyOrbit{T<:AbstractConic} <: TwoBodySystem{T}

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
    InvalidTwoBodyOrbit(body::CelestialBody)

Returns a `TwoBodyOrbit` with `NaN` state values. Used by 
`propagate_twobody` and `kepler` to indicate failed convergance.
"""
InvalidTwoBodyOrbit(body::CelestialBody) = TwoBodyOrbit(SVector{3}(NaN * u"km", NaN * u"km", NaN * u"km"),
                                                        SVector{3}(NaN * u"km", NaN * u"km", NaN * u"km"),
                                                        body)

"""
    isinvalid(orbit::TwoBodyOrbit)

Checks for `NaN` valued orbital states, which are used to
indicate an invalid `TwoBodyOrbit`.
"""
isinvalid(orbit::TwoBodyOrbit)    = all(map(x->!isnan(getfield(orbit, x)), [:r̅, :v̅, :e, :a, :i, :Ω, :ω, :ν])) 