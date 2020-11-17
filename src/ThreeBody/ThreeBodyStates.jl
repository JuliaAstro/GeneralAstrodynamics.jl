#
# Handles CR3BP problem states.
#

"""
    ThreeBodySystem

Describes a Circular Restricted Three-Body
Problem system.
"""
struct ThreeBodySystem{F<:AbstractFloat} <: OrbitalSystem

    # Dimensional Units
    μ₁::Unitful.Quantity{F}
    μ₂::Unitful.Quantity{F}
    r::Unitful.Length{F}
    v::Unitful.Velocity{F}
    t::Unitful.Time{F}

    # Non-dimensional Units
    M::Unitful.DimensionlessQuantity{F}
    μ::Unitful.DimensionlessQuantity{F}
    x₁::Unitful.DimensionlessQuantity{F}
    x₂::Unitful.DimensionlessQuantity{F}
    T::Unitful.DimensionlessQuantity{F}
    
end