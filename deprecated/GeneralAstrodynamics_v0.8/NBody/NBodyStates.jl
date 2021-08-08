#
# Handles NBody problem states.
#

"""
Stores the state of each body in the NBody problem.
"""
struct Body{F<:AbstractFloat} <: AbstractBody
   
    r::SVector{3, Unitful.Length{F}}
    v::SVector{3, Unitful.Velocity{F}}
    m::Unitful.Mass{F}
    name::String

    function Body(r::VR, v::VV, m::M, name::String="") where {
            T1  <: Real, 
            T2  <: Real, 
            T3  <: Real, 
            R   <: Unitful.Length{T1},
            V   <: Unitful.Velocity{T2},
            VR  <: AbstractVector{R}, 
            VV  <: AbstractVector{V},
            M   <: Unitful.Mass{T3}
        }
        if length(r) ≢ length(v) ≢ 3
            throw(ArgumentError("The `Body` constructor requires 3 element vectors for position `r` and velocity `v`"))
        else
            T = promote_type(T1, T2, T3)
            if !(T <: AbstractFloat)
                @warn "Non-float parameters provided. Defaulting to Float64."
                T = Float64
            end
            return new{T}(SVector{3}(T.(r)), SVector{3}(T.(v)), T(m), name)
        end
    end

    function Body(r::VR, v::VV, m::M, name::String="") where {
            R   <: Real,
            V   <: Real,
            VR  <: AbstractVector{R}, 
            VV  <: AbstractVector{V},
            M   <: Real
        }
        @warn "No units provided! Assuming km, km/s, and kg."
        return Body(r * u"km", v * u"km/s", m * u"kg", name)
    end
end

Base.convert(::Type{T}, b::Body) where {T<:AbstractFloat} = Body(T.(b.r), T.(b.v), T(b.m), b.name)
Base.promote(::Type{Body{A}}, ::Type{Body{B}}) where {A<:AbstractFloat, B<:AbstractFloat} = Body{promote_type(A,B)}
Core.Float16(o::Body) = convert(Float16, o)
Core.Float32(o::Body) = convert(Float32, o)
Core.Float64(o::Body) = convert(Float64, o)
Base.MPFR.BigFloat(o::Body) = convert(BigFloat, o)

"""
Describes a system of `n` `NBodyStates`'s.
"""
struct NBodySystem{N, T<:AbstractFloat} <: AbstractOrbitalSystem

    bodies::SVector{N, Body{T}}

    function NBodySystem(b::B...) where B<:Body
        N = length(b)
        bodies = promote(b...)
        T = promote_type([typeof(bᵢ).parameters[1] for bᵢ ∈ b]...)
        return new{length(b),T}(SVector{N, Body{T}}(bodies...))
    end
    NBodySystem(b::VB) where {B<:Body, VB<:AbstractVector{B}} = NBodySystem(b...)

end

"""
The `length` of an `NBodySystem` is the number of bodies in the system.
"""
Base.@pure Base.length(sys::NBodySystem{N,T}) where N where T = N

"""
The n-th `index` of an `NBodySystem` is the n-th body in the system.
"""
Base.getindex(sys::NBodySystem, i) = sys.bodies[i]
Base.convert(::Type{T}, sys::NBodySystem) where {T<:AbstractFloat} = NBodySystem(convert.(T, sys.bodies)...)
Base.promote(::Type{NBodySystem{N, A}}, ::Type{NBodySystem{N, B}}) where {A<:AbstractFloat, B<:AbstractFloat,N} = NBodySystem{promote_type(A,B)}
Core.Float16(o::NBodySystem) = convert(Float16, o)
Core.Float32(o::NBodySystem) = convert(Float32, o)
Core.Float64(o::NBodySystem) = convert(Float64, o)
Base.MPFR.BigFloat(o::NBodySystem) = convert(BigFloat, o)
