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
    m₁::Unitful.Mass{F}
    m₂::Unitful.Mass{F}
    x₁::Unitful.Length{F}
    x₂::Unitful.Length{F}

    # Non-dimensional Units
    μ::Unitful.Quantity{F}
    
    
end