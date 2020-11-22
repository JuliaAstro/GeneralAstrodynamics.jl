#
# Handles CR3BP problem states.
#

"""
Describes a Circular Restricted Three-Body
Problem system.
"""
struct ThreeBodySystem{F<:AbstractFloat} <: OrbitalSystem

    # Dimensional Units
    a::Length{F}
    μ₁::MassParameter{F}
    μ₂::MassParameter{F}
    t::Time{F}

    # Non-dimensional Units
    μ::F
    rₛ::SVector{3, F}
    vₛ::SVector{3, F}
    tₛ::F

    function ThreeBodySystem(a::A, μ₁::U1, μ₂::U2, 
                             r::VR, v::VV,  t::TT) where {
            A  <: Length,
            U1 <: MassParameter,
            U2 <: MassParameter,
            R  <: Length, 
            V  <: Velocity, 
            VR <: AbstractVector{R}, 
            VV <: AbstractVector{V},
            TT <: Time
        }

        T = typeof(promote(a, μ₁, μ₂, r, v, t)[1]).val

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
        return new{T}(map(x->T.(x), (
            a, μ₁, μ₂, SVector{3}(r), SVector{3}(v), t, Tₛ, 
            μ, -μ, 1-μ, nondimensionalize(SVector{3}(r), a), nondimensionalize(SVector{3}(v), a, Tₛ), t/Tₛ))...
        )

    end

end