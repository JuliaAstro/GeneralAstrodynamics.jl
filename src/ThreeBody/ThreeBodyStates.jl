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

function Base.show(io::IO, sys::ThreeBodyState)
    println(io, "Dimensioned Circular Restricted Three-body State")
    println(io, "  μ₁:        ", sys.μ₁)
    println(io, "  μ₂:        ", sys.μ₂)
    println(io, "   a:        ", sys.a)
    println(io, "  r₃:        ", transpose(ustrip.(sys.r₃)), " ", unit(first(sys.r₃)))
    println(io, "  v₃:        ", transpose(ustrip.(sys.v₃)), " ", unit(first(sys.v₃)))
    println(io, "  Δt:        ", sys.Δt)
end

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

    function NondimensionalThreeBodyState(rₛ::AbstractVecOrMat{R}, vₛ::AbstractVecOrMat{V}, μ::U, Δt::D = 1.0, 
                                          DU::Unitful.Length{L} = NaN * u"km", 
                                          DT::Unitful.Time{C} = NaN * u"s") where {
                                          R <: Real, V <: Real, L <: Real, C <: Real, U <: Real, D <: Real}
        T = promote_type(R, V, L, C, U, D)
        if !(T <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            T = Float64
        end
        
        return new{T}(
            SVector{3,T}(rₛ...), 
            SVector{3,T}(vₛ...), 
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

function Base.show(io::IO, sys::NondimensionalThreeBodyState)
    println(io, "Nondimensional Circular Restricted Three-body State")
    println(io, "   μ:        ", sys.μ)
    println(io, "   r:        ", transpose(sys.r))
    println(io, "   v:        ", transpose(sys.v))
    println(io, "  Δt:        ", sys.Δt)
    println(io, "  DU:        ", sys.DU)
    println(io, "  DT:        ", sys.DT)
end

"""
Returns true if all elements in each system are within `atol` of the other.
"""
function Base.isapprox(c1::ThreeBodyState, c2::ThreeBodyState; atol = 1e-8)
    ru(x) = ustrip(upreferred(x))
    return isapprox(ru(c1.μ₁),  ru(c2.μ₁);  atol = atol) &&
           isapprox(ru(c1.μ₂),  ru(c2.μ₂);  atol = atol) &&
           isapprox(ru(c1.a),   ru(c2.a);   atol = atol) &&
           isapprox(ru.(c1.r₃), ru.(c2.r₃); atol = atol) &&
           isapprox(ru.(c1.v₃), ru.(c2.v₃); atol = atol) &&
           isapprox(ru(c1.Δt),  ru(c2.Δt);  atol = atol)

end

"""
Returns true if all elements in each system are equal to the other.
"""
function Base.isequal(c1::ThreeBodyState, c2::ThreeBodyState)
    ru(x) = ustrip(upreferred(x))
    return isequal(ru(c1.μ₁),  ru(c2.μ₁))  &&
           isequal(ru(c1.μ₂),  ru(c2.μ₂))  &&
           isequal(ru(c1.a),   ru(c2.a))   &&
           isequal(ru.(c1.r₃), ru.(c2.r₃)) &&
           isequal(ru.(c1.v₃), ru.(c2.v₃)) &&
           isequal(ru(c1.Δt),  ru(c2.Δt))

end

"""
Returns true if all elements in each system are within `atol` of the other.
"""
function Base.isapprox(c1::NondimensionalThreeBodyState, c2::NondimensionalThreeBodyState; atol = 1e-8)
    ru(x) = ustrip(upreferred(x))
    return isapprox(ru(c1.μ),   ru(c2.μ);   atol = atol) &&
           isapprox(ru.(c1.r),  ru.(c2.r);  atol = atol) &&
           isapprox(ru.(c1.v),  ru.(c2.v);  atol = atol) &&
           isapprox(ru(c1.Δt),  ru(c2.Δt);  atol = atol) &&
           isapprox(ru.(c1.DU), ru(c2.DU);  atol = atol) &&
           isapprox(ru(c1.DT),  ru(c2.DT);  atol = atol)

end

"""
Returns true if all elements in each system are equal to the other.
"""
function Base.isequal(c1::NondimensionalThreeBodyState, c2::NondimensionalThreeBodyState)
    ru(x) = ustrip(upreferred(x))
    return isequal(ru(c1.μ),  ru(c2.μ))  &&
           isequal(ru.(c1.r), ru.(c2.r)) &&
           isequal(ru.(c1.v), ru.(c2.v)) &&
           isequal(ru(c1.Δt), ru(c2.Δt)) &&
           isequal(ru(c1.DU), ru(c2.DU)) &&
           isequal(ru(c1.DT),  ru(c2.DT))

end