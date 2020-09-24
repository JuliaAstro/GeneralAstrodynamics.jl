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
abstract type TwoBodyState{T<:AbstractConic} <: AbstractOrbit end


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

    CartesianState(r̅, v̅, body)

Cartesian representation for orbital state.
"""
struct CartesianState{T<:AbstractConic} <: TwoBodyState{T}

    r̅::SVector{3, Unitful.Length}
    v̅::SVector{3, Unitful.Velocity}
    body::Body

end

"""
    KeplerianState(e, a, i , Ω, ω, ν, body)

Keplarian representation for orbital state.
"""
struct KeplerianState{T<:AbstractConic} <: TwoBodyState{T}

    e::Unitful.DimensionlessQuantity
    a::Unitful.Length
    i::Unitful.DimensionlessQuantity
    Ω::Unitful.DimensionlessQuantity
    ω::Unitful.DimensionlessQuantity
    ν::Unitful.DimensionlessQuantity
    body::Body

end

"""
    TwoBodyOrbit{T<:AbstractConic}

Struct for storing TwoBody orbital states for all conics.
"""
struct TwoBodyOrbit{T<:AbstractConic} <: TwoBodyState{T}

    # Cartesian representation
    r̅::SVector{3, Unitful.Length}
    v̅::SVector{3, Unitful.Velocity}

    # Keplerian representation
    e::Unitful.DimensionlessQuantity
    a::Unitful.Length
    i::Unitful.DimensionlessQuantity
    Ω::Unitful.DimensionlessQuantity
    ω::Unitful.DimensionlessQuantity
    ν::Unitful.DimensionlessQuantity

    # Body
    body::Body

end