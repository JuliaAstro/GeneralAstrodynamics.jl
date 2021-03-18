#
# Handles CR3BP problem states.
#

"""
Abstract type for restricted three-body systems.
"""
abstract type RestrictedThreeBodySystem <: OrbitalSystem end

"""
Describes a dimensional state of a spacecraft
within the Circular Restriested Three-body Problem in
the Synodic frame.
"""
struct ThreeBodyState{F<:AbstractFloat} <: RestrictedThreeBodySystem

    μ₁::MassParameter{F}
    μ₂::MassParameter{F}
    a::Length{F}
    r₃::SVector{3, <:Length{F}}
    v₃::SVector{3, <:Velocity{F}}
    Δt::Time{F}

    function ThreeBodyState(μ₁::MP1, μ₂::MP2, a::A, r₃::R, v₃::V, Δt::DT) where {
            MP1 <: MassParameter{<:Real},
            MP2 <: MassParameter{<:Real},
            A   <: Length{<:Real},
            RT  <: Length{<:Real},
            VT  <: Velocity{<:Real},
            R   <: AbstractVector{RT},
            V   <: AbstractVector{VT},
            DT  <: Time{<:Real}
        }

        T = promote_type(
            typeof(μ₁.val), typeof(μ₂.val), typeof(a.val),
            map(x->typeof(x.val), r₃)..., map(x->typeof(x.val), v₃)...,
            typeof(Δt.val)
        )
        if !(T <: AbstractFloat)
            @warn "Non-float parameters provided: defaulting to Float64."
            T = Float64
        end

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
Core.Float16(o::ThreeBodyState) = convert(Float16, o)
Core.Float32(o::ThreeBodyState) = convert(Float32, o)
Core.Float64(o::ThreeBodyState) = convert(Float64, o)
Base.MPFR.BigFloat(o::ThreeBodyState) = convert(BigFloat, o)

"""
Describes the non-dimensional state of a spacecraft
within the Circular Restricted Three-body Problem in
the Synodic frame.
"""
struct NondimensionalThreeBodyState{F<:AbstractFloat} <: RestrictedThreeBodySystem
    
    r::SVector{3, F}
    v::SVector{3, F}
    μ::F
    Δt::F
    DU::Length{F}
    DT::Time{F}

    function NondimensionalThreeBodyState(rₛ::R, vₛ::V, μ::F1, Δt::F2 = 1.0, 
                                          DU::Unitful.Length{F3} = NaN * u"km", 
                                          DT::Unitful.Time{F4} = NaN * u"km/s") where {
            FR <: Real,
            FV <: Real,
            R  <: AbstractArray{FR}, 
            V  <: AbstractArray{FV},
            F1 <: Real,
            F2 <: Real,
            F3 <: Real,
            F4 <: Real
        }

        T = promote_type(FR, FV, F1, F2, F3, F4)
        if !(T <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            T = Float64
        end
        
        return new{T}(
            SVector{3,T}(T.(rₛ[:])), 
            SVector{3,T}(T.(vₛ[:])), 
            T(μ), T(Δt), T(DU), T(DT)
        )

    end

end

Base.convert(::Type{T}, t::NondimensionalThreeBodyState) where {
    T<:AbstractFloat
} = NondimensionalThreeBodyState(T.(Array(t.rₛ)), T.(Array(t.vₛ)), T(t.μ), T(t.Δt), T(t.DU), T(t.DT))
Base.promote(::Type{NondimensionalThreeBodyState{A}}, ::Type{NondimensionalThreeBodyState{B}}) where {
    A<:AbstractFloat, B<:AbstractFloat
} = NondimensionalThreeBodyState{promote_type(A,B)}
Core.Float16(o::NondimensionalThreeBodyState) = convert(Float16, o)
Core.Float32(o::NondimensionalThreeBodyState) = convert(Float32, o)
Core.Float64(o::NondimensionalThreeBodyState) = convert(Float64, o)
Base.MPFR.BigFloat(o::NondimensionalThreeBodyState) = convert(BigFloat, o)