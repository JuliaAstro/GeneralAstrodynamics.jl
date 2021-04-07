"""
Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module CommonTypes

using Reexport
@reexport using Unitful, UnitfulAngles, UnitfulAstro
using StaticArrays: StaticVector, MVector

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

export AbstractBody, AbstractOrbitalSystem, AbstractTrajectory, AbstractCartesianState, CartesianState
export getindex, setindex!, lengthunit, timeunit, velocityunit
export position_vector, velocity_vector, scalar_position, scalar_velocity
export convert, promote, Float16, Float32, Float64, BigFloat

""" 
Abstract type for bodies in space: both `CelestialBody`s (in
`TwoBody.jl`), and `Body`s (in `NBody.jl`).
"""
abstract type AbstractBody end

"""
Abstract type describing select astrodynamics problems.
"""
abstract type AbstractOrbitalSystem end

"""
Abstract type describing a collection of states resulting from numerical integration
"""
abstract type AbstractTrajectory end

"""
Abstract type describing an orbital state.
"""
abstract type AbstractCartesianState{F<:AbstractFloat} <: StaticVector{6,F} end

"""
Cartesian state which describes a spacecraft or body's position and velocity with respect to _something_.
"""
mutable struct CartesianState{F<:AbstractFloat, LU, TU} <: AbstractCartesianState{F} where {LU<:Unitful.LengthUnits, TU<:Unitful.TimeUnits}
    r::SubArray{F, 1, MVector{6, F}, Tuple{UnitRange{Int64}}, true}
    v::SubArray{F, 1, MVector{6, F}, Tuple{UnitRange{Int64}}, true}
    rv::MVector{6,F}

    function CartesianState(r::AbstractVector{R}, v::AbstractVector{V}; lengthunit=u"km", timeunit=u"s") where {R<:Real, V<:Real}
        F  = promote_type(R, V)
        if !(F <: AbstractFloat)
            @warn "Type provided ($(string(F))) is not a float: defaulting to Float64."
            F = Float64
        end
        @assert length(r) == length(v) == 3 "Both arguments must have length equal to 3!"
        rv = MVector{6, F}(r[1], r[2], r[3], v[1], v[2], v[3])
        pos = @views rv[1:3]
        vel = @views rv[4:6]
        return new{F, lengthunit, timeunit}(pos, vel, rv)
    end

    function CartesianState(r::AbstractVector{R}, v::AbstractVector{V}) where {R<:Unitful.Length, V<:Unitful.Velocity}

        # Get ready to commit a crime... we need the Time unit provided in velocity quantity V

        # V is a quantity with units Length / Time.  Let's reverse that so we have Time instead of Time^(-1)
        TL = inv(unit(V))

        # Now we need to find which index (1 or 2) contains the Time unit
        timeaxis = findfirst(T -> T isa Unitful.Dimension{:Time}, collect(typeof(typeof(TL).parameters[2]).parameters[1]))
        
        # Now that we have the proper index, let's select the time unit
        timeunit = typeof(TL).parameters[1][timeaxis]

        # Phew!
        return CartesianState(ustrip.(unit(R), r), ustrip.(unit(V), v); lengthunit = unit(R), timeunit = timeunit)
    end
end

Base.convert(::Type{T}, o::CartesianState) where {T<:AbstractFloat} = CartesianState(T.(o.r), T.(o.v))
Base.promote(::Type{CartesianState{A}}, ::Type{CartesianState{B}}) where {A<:AbstractFloat, B<:AbstractFloat} = CartesianState{promote_type(A,B)}
Core.Float16(o::CartesianState) = convert(Float16, o)
Core.Float32(o::CartesianState) = convert(Float32, o)
Core.Float64(o::CartesianState) = convert(Float64, o)
Base.MPFR.BigFloat(o::CartesianState) = convert(BigFloat, o)

Base.getindex(state::CartesianState, i::Int) = state.rv[i]
Base.setindex!(state::CartesianState, value, i::Int) = (state.rv[i] = value)

"""
Returns the `Unitful.Length` unit associated with the Cartesian state.
"""
Base.@pure lengthunit(::C) where C <: CartesianState = C.parameters[2]

"""
Returns the `Unitful.Time` unit associated with the Cartesian state.
"""
Base.@pure timeunit(::C) where C <: CartesianState = C.parameters[3]

"""
Returns the `Unitful.Velocity` unit associated with the Cartesian state.
"""
velocityunit(state::CartesianState) = lengthunit(state) / timeunit(state)

function Base.show(io::IO, ::MIME"text/plain", X::CartesianState{F, LU, TU}) where {F, LU, TU} 
    println(io, "Cartesian State of type $(string(F)):")
    println(io, "  r = ", [state.r[1] state.r[2] state.r[3]], " ", string(LU))
    println(io, "  v = ", [state.v[1] state.v[2] state.v[3]], " ", string(LU/TU))
end

"""
Returns the `Unitful` position vector of the `Cartesianstate`.
"""
position_vector(state::CartesianState) = state.r * lengthunit(state)

"""
Returns the `Unitful` velocity vector of the `CartesianState`.
"""
velocity_vector(state::CartesianState) = state.v * velocityunit(state)

"""
Returns the `Unitful` scalar position of the `CartesianState`.
"""
scalar_position(state::CartesianState) = norm(position_vector(state))

"""
Returns the `Unitful` scalar velocity of the `CartesianState`.
"""
scalar_velocity(state::CartesianState) = norm(velocity_vector(state))

end