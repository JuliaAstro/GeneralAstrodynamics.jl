#
# Handles NBody problem states.
#

"""
Stores the state of each body in the NBody problem.
"""
struct Body{F<:AbstractFloat} <: AbstractBody
   
    r̅::SVector{3, Unitful.Length{F}}
    v̅::SVector{3, Unitful.Velocity{F}}
    m::Unitful.Mass{F}

    function Body(r::VR, v::VV, m::M) where {
            T1  <: Real, 
            T2  <: Real, 
            T3  <: Real, 
            R   <: Unitful.Length{T1},
            V   <: Unitful.Velocity{T2},
            VR  <: AbstractVector{R}, 
            VV  <: AbstractVector{V},
            M   <: Unitful.Mass{T3}}

        if length(r) ≢ length(v) ≢ 3
            error("The `Body` constructor requires 3 element vectors for position `r` and velocity `v`")
        else
            T = promote_type(T1, T2, T3)
            if !(T <: AbstractFloat)
                @warn "Non-float types provided to `Body` constructor: using `Float64`."
                T = Float64
            end
            return new{T}(SVector{3}(T.(r)), SVector{3}(T.(v)), T(m))
        end

    end

end
Base.convert(::Type{T}, b::Body) where {
        T<:AbstractFloat} = Body(T.(b.r̅), T.(b.v̅), T(b.m))
Base.promote(::Type{Body{A}}, ::Type{Body{B}}) where {
        A<:AbstractFloat, B<:AbstractFloat} = Body{promote_type(A,B)}



"""
Describes a system of `n` `MultibodyStates`'s.
"""
struct MultibodySystem{N, T<:AbstractFloat} <: OrbitalSystem

    body::SVector{N, Body{T}}

    function MultibodySystem(b::B...) where B<:Body
        N = length(b)
        bodies = promote(b...)
        T = typeof(first(bodies).m.val)
        return new{length(b),T}(SVector{N, Body{T}}(bodies...))
    end
    MultibodySystem(b::VB) where {B<:Body, VB<:AbstractVector{B}} = MultibodySystem(b...)

end
Base.convert(::Type{T}, sys::MultibodySystem) where {
        T<:AbstractFloat} = MultibodySystem(convert.(Ref(T), sys.body)...)
Base.promote(::Type{MultibodySystem{N, A}}, ::Type{MultibodySystem{N, B}}) where {
        A<:AbstractFloat, B<:AbstractFloat,N} = MultibodySystem{promote_type(A,B)}

