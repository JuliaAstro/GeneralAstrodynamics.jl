#
# Handles CR3BP problem states.
#

"""
Describes a dimensional state of a spacecraft
within the Circular Restriested Three-body Problem in
the Synodic frame.
"""
struct ThreeBodyState{F<:AbstractFloat} <: OrbitalSystem

    μ₁::MassParameter{F}
    μ₂::MassParameter{F}
    a::Length{F}
    r₃::SVector{3, <:Length{F}}
    v₃::SVector{3, <:Velocity{F}}
    Δt::Time{F}

    function ThreeBodyState(μ₁::MP1, μ₂::MP2, a::A, r₃::R, v₃::V, Δt::DT) where {
            MP1 <: MassParameter{<:AbstractFloat},
            MP2 <: MassParameter{<:AbstractFloat},
            A   <: Length{<:AbstractFloat},
            RT  <: Length{<:AbstractFloat},
            VT  <: Velocity{<:AbstractFloat},
            R   <: AbstractVector{RT},
            V   <: AbstractVector{VT},
            DT  <: Time{<:AbstractFloat}
        }

        T = promote_type(
            typeof(μ₁.val), typeof(μ₂.val), typeof(a.val),
            map(x->typeof(x.val), r₃)..., map(x->typeof(x.val), v₃)...,
            typeof(Δt.val)
        )

        return new{T}(
            T(μ₁), T(μ₂), T(a), 
            SVector{3, RT}(T.(r₃)),
            SVector{3, VT}(T.(v₃)),
            T(Δt)
        )

    end

end
Base.convert(::Type{T}, t::ThreeBodyState) where {
        T<:AbstractFloat
    } = ThreeBodyState(map(f->T.(getfield(t,f), fieldnames(t)))...)
Base.promote(::Type{ThreeBodyState{A}}, ::Type{ThreeBodyState{B}}) where {
        A<:AbstractFloat, B<:AbstractFloat
    } = ThreeBodyState{promote_type(A,B)}

"""
Describes the non-dimensional state of a spacecraft
within the Circular Restricted Three-body Problem in
the Synodic frame.
"""
struct NondimensionalThreeBodyState{F<:AbstractFloat} <: OrbitalSystem
    
    r₁::SVector{3, F}
    r₂::SVector{3, F}
    rₛ::SVector{3, F}
    vₛ::SVector{3, F}
    μ::F
    DU::Length{F}
    DT::Time{F}

    function NondimensionalThreeBodyState(rₛ::R, vₛ::V, μ::F1, DU::F2, DT::F3) where {
            R  <: AbstractArray{<:AbstractFloat}, 
            V  <: AbstractArray{<:AbstractFloat},
            F1 <: AbstractFloat,
            F2 <: Length{<:AbstractFloat},
            F3 <: Time{<:AbstractFloat}
        }

        T = promote_type(eltype(rₛ), eltype(vₛ), typeof(μ), typeof(DU.val), typeof(DT.val))
        return new{T}(
            SVector{3,T}(T.([-μ,  0, 0])),
            SVector{3,T}(T.([1-μ, 0, 0])),
            SVector{3,T}(T.(rₛ[:])), 
            SVector{3,T}(T.(vₛ[:])), 
            T(μ), T(DU), T(DT)
        )

    end

end
Base.convert(::Type{T}, t::NondimensionalThreeBodyState) where {
        T<:AbstractFloat
    } = NondimensionalThreeBodyState(T.(Array(t.rₛ)), T.(Array(t.vₛ)), T(t.μ), T(t.DU), T(t.DT))
Base.promote(::Type{NondimensionalThreeBodyState{A}}, ::Type{NondimensionalThreeBodyState{B}}) where {
        A<:AbstractFloat, B<:AbstractFloat
    } = NondimensionalThreeBodyState{promote_type(A,B)}

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

        T = promote_type(typeof(a.val), 
                         typeof(μ₁.val), 
                         typeof(μ₂.val), 
                         map(x->typeof(x.val), r)...,
                         map(x->typeof(x.val), v)...,
                         typeof(t.val))

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
