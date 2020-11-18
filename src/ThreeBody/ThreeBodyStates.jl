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
    a::Unitful.Length{F}
    μ₁::Unitful.Quantity{F}
    μ₂::Unitful.Quantity{F}
    r::SVector{3, Unitful.Length{F}}
    v::SVector{3, Unitful.Velocity{F}}
    t::Unitful.Time{F}
    T::Unitful.Time{F}

    # Non-dimensional Units
    μ::F
    x₁::F
    x₂::F
    rₛ::SVector{3, F}
    vₛ::SVector{3, F}
    tₛ::F

    function ThreeBodySystem(a, μ₁, μ₂, r, v, t)

        if length(r) ≢ length(v) ≢ 3
        throw(ArgumentError(string("Both `r` and `v` provided to `ThreeBodySystem` ",
            "constructor must have length 3.")))
        end

        if μ₂ > μ₁
        @warn string("The second mass parameter is larger than the first. ",
                     "Did you mean to switch those two? Assuming the second ",
                     "mass parameter is the primary body.")
        end

        μ = min(μ₁, μ₂) / (μ₁+μ₂)
        Tₛ = orbital_period(a, μ₁+μ₂)
        return new{Float64}(map(x->Float64.(x), (
            a, μ₁, μ₂, SVector{3}(r), SVector{3}(v), t, Tₛ, 
            μ, -μ, 1-μ, nondimensionalize(SVector{3}(r), a), nondimensionalize(SVector{3}(v), a, Tₛ), t/Tₛ))...
        )

    end

    function ThreeBodySystem( a::Unitful.Length{T}, 
                             μ₁::Unitful.Quantity{T},
                             μ₂::Unitful.Quantity{T}, 
                              r::VR, 
                              v::VV, 
                              t::Unitful.Time{T}) where {T  <: AbstractFloat, 
                                                         R  <: Unitful.Length{T}, 
                                                         V  <: Unitful.Velocity{T}, 
                                                         VR <: AbstractVector{R}, 
                                                         VV <: AbstractVector{V}}

        if length(r) ≢ length(v) ≢ 3
            throw(ArgumentError(string("Both `r` and `v` provided to `ThreeBodySystem` ",
                                       "constructor must have length 3.")))
        end

        if μ₂ > μ₁
            @warn string("The second mass parameter is larger than the first. ",
                         "Did you mean to switch those two? Assuming the second ",
                         "mass parameter is the primary body.")
        end

        μ = min(μ₁, μ₂) / (μ₁+μ₂)
        Tₛ = orbital_period(a, μ₁+μ₂)
        return new{T}(
            a, μ₁, μ₂, SVector{3}(r), SVector{3}(v), t, Tₛ, 
            μ, -μ, 1-μ, nondimensionalize(SVector{3}(r), a), nondimensionalize(SVector{3}(v), a, Tₛ), nondimensionalize(t, Tₛ)
        )

    end

end