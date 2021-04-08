"""
Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module AstrodynamicsCore

export  AbstractBody, AbstractOrbitalState, AbstractUnitfulState, AbstractOrbitalSystem, AbstractTrajectory, AbstractUnitfulState, CartesianState,
        getindex, setindex!, lengthunit, timeunit, velocityunit, massparameterunit,
        position_vector, velocity_vector, scalar_position, scalar_velocity

using Reexport
@reexport using Unitful, UnitfulAngles, UnitfulAstro
using StaticArrays: StaticVector, MVector

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

""" 
Abstract type for bodies in space: both `CelestialBody`s (in
`TwoBody.jl`), and `Body`s (in `NBody.jl`).
"""
abstract type AbstractBody end

"""
Abstract type describing select astrodynamics problems.
"""
abstract type AbstractOrbitalState end

"""
Abstract type describing select astrodynamics problem constants.
"""
abstract type AbstractOrbitalSystem end

"""
Abstract type describing a collection of states resulting from numerical integration
"""
abstract type AbstractTrajectory end

"""
Abstract type describing an orbital state.
"""
abstract type AbstractUnitfulState{F<:AbstractFloat, LU, TU} <: StaticVector{6,F} where {LU <: Unitful.LengthFreeUnits, TU <: Unitful.TimeFreeUnits} end

"""
Cartesian state which describes a spacecraft or body's position and velocity with respect to _something_.
"""
mutable struct CartesianState{F, LU, TU} <: AbstractUnitfulState{F, LU, TU}
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
        timeunit = Unitful.Units{(typeof(TL).parameters[1][timeaxis],), Unitful.𝐓}()

        # This is easy...
        lengthunit = unit(R)

        # Phew!
        return CartesianState(ustrip.(lengthunit, r), ustrip.(lengthunit/timeunit, v); lengthunit = lengthunit, timeunit = timeunit)
    end

    function CartesianState(cart::AbstractUnitfulState{F}) where {F}
        return CartesianState(cart.r,  cart.v; lengthunit = lengthunit(cart), timeunit = timeunit(cart))
    end

    function CartesianState(arr::StaticVector{6})
        return CartesianState(arr[1:3], arr[4:6])
    end
end

Base.getindex(cart::CartesianState, i::Int) = cart.rv[i]
Base.setindex!(cart::CartesianState, value, i::Int) = (cart.rv[i] = value)

"""
Returns the `Unitful.Length` unit associated with the state.
"""
lengthunit(::S) where S <: AbstractUnitfulState = S.parameters[2]

"""
Returns the `Unitful.Time` unit associated with the state.
"""
timeunit(::S) where S <: AbstractUnitfulState = S.parameters[3]

"""
Returns the `Unitful.Velocity` unit associated with the state.
"""
velocityunit(state::S) where S <: AbstractUnitfulState = S.parameters[2] / S.parameters[3]

"""
Returns the `MassParameter` unit associated with the state.
"""
massparameterunit(state::S) where S <: AbstractUnitfulState = S.parameters[2]^3 / S.parameters[3]^2

function Base.show(io::IO, state::CartesianState{F, LU, TU}) where {F<:AbstractFloat, LU, TU} 
    println(io, "  Cartesian State:")
    println(io, "    r = ", [state.r[1] state.r[2] state.r[3]], " ", string(LU))
    println(io, "    v = ", [state.v[1] state.v[2] state.v[3]], " ", string(LU/TU))
end

function Base.show(io::IO, ::MIME"text/plain", state::CartesianState{F, LU, TU}) where {F<:AbstractFloat, LU, TU} 
    println(io, "  Cartesian State:")
    println(io, "    r = ", [state.r[1] state.r[2] state.r[3]], " ", string(LU))
    println(io, "    v = ", [state.v[1] state.v[2] state.v[3]], " ", string(LU/TU))
end

"""
Returns the `Unitful` position vector of the `Cartesianstate`.
"""
position_vector(state::CartesianState{F,LU,TU}) where {F, LU, TU} = state.r * lengthunit(state)

"""
Returns the `Unitful` velocity vector of the `CartesianState`.
"""
velocity_vector(state::CartesianState{F,LU,TU}) where {F, LU, TU} = state.v * velocityunit(state)

"""
Returns the `Unitful` scalar position of the `CartesianState`.
"""
scalar_position(state::CartesianState{F,LU,TU}) where {F, LU, TU} = norm(position_vector(state))

"""
Returns the `Unitful` scalar velocity of the `CartesianState`.
"""
scalar_velocity(state::CartesianState{F,LU,TU}) where {F, LU, TU} = norm(velocity_vector(state))

end