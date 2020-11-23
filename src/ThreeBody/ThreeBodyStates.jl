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
    rₛ::SVector{3, F}
    vₛ::SVector{3, F}
    tₛ::F
    μ::F

    function ThreeBodySystem(a::A, μ₁::U1, μ₂::U2, 
                             r::VR, v::VV,  t::TT) where {
            A  <: Length{<:AbstractFloat},
            U1 <: MassParameter{<:AbstractFloat},
            U2 <: MassParameter{<:AbstractFloat},
            R  <: Length{<:AbstractFloat}, 
            V  <: Velocity{<:AbstractFloat}, 
            VR <: AbstractVector{R}, 
            VV <: AbstractVector{V},
            TT <: Time{<:AbstractFloat}
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

        return new{T}(a, μ₁, μ₂, t, nondimensionalize(r, v, t, μ₁, μ₂, a)...)

    end
    ThreeBodySystem(a, μ₁, μ₂, r, v, t) = ThreeBodySystem(
        Float64(a), Float64(μ₁), Float64(μ₂),
        Float64.(r), Float64.(v), Float64(t))

end
Base.convert(::Type{T}, t::ThreeBodySystem) where {
        T<:AbstractFloat
    } = ThreeBodySystem(map(f->T.(getfield(t,f), fieldnames(t)))...)
Base.promote(::Type{ThreeBodySystem{A}}, ::Type{ThreeBodySystem{B}}) where {
        A<:AbstractFloat, B<:AbstractFloat
    } = ThreeBodySystem{promote_type(A,B)}
