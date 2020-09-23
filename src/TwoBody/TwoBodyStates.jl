#
#   TwoBodyStates.jl
#
#   Describes Two Body Orbits through Cartesian coordinates and Orbital Elements.
# 

"""
    TwoBodyOrbit

Abstract type for all two-body orbital representations.
"""
abstract type TwoBodyOrbit <: AbstractOrbit end

"""
    CartesianState(r̅, v̅, body)

Cartesian representation for orbital state.
"""
struct CartesianState <: TwoBodyOrbit

    r̅::SVector{3, Unitful.Length}
    v̅::SVector{3, Unitful.Velocity}
    body::Body

end

"""
    KeplerianState(e, a, i , Ω, ω, ν, body)

Keplarian representation for orbital state.
"""
struct KeplerianState <: TwoBodyOrbit

    e::Unitful.DimensionlessQuantity
    a::Unitful.Length
    i::Unitful.Quantity
    Ω::Unitful.Quantity
    ω::Unitful.Quantity
    ν::Unitful.Quantity
    body::Body

end

