"""
Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module AstrodynamicsCore

export  AbstractBody, AbstractOrbitalSystem, AbstractUnitfulStructure, AbstractState, AbstractSystem, AbstractOrbit, AbstractTrajectory, CartesianState,
        getindex, setindex!, lengthunit, timeunit, velocityunit, massparameterunit, coordinateframe,
        position_vector, velocity_vector, scalar_position, scalar_velocity,
        AbstractFrame, Bodycentric, Synodic, Perifocal, size, length, convert, epoch, eltype

export NormalizedLengthUnit, NormalizedTimeUnit

using Reexport
@reexport using Unitful, UnitfulAngles, UnitfulAstro
using StaticArrays: StaticVector, MVector
using LinearAlgebra: norm

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

""" 
Abstract type for bodies in space: both `CelestialBody`s (in
`TwoBody.jl`), and `Body`s (in `NBody.jl`).
"""
abstract type AbstractBody end


abstract type AbstractOrbitalSystem end

"""
Abstract type describing a collection of states resulting from numerical integration
"""
abstract type AbstractTrajectory end

"""
Absract type describing all frames.
"""
abstract type AbstractFrame end

abstract type AbstractUnitfulStructure{F <: AbstractFloat, LU <: Unitful.LengthUnits, TU <: Unitful.TimeUnits} end

abstract type AbstractState{F, LU, TU, FR} <: AbstractUnitfulStructure{F, LU, TU} where FR <: AbstractFrame end

abstract type AbstractSystem{F, LU, TU} <: AbstractUnitfulStructure{F, LU, TU} end

abstract type AbstractOrbit{F, LU, TU} <: AbstractUnitfulStructure{F, LU, TU} end

"""
Body centered frame.
"""
struct Bodycentric <: AbstractFrame end

"""
Synodic frame.
"""
struct Synodic <: AbstractFrame end

"""
Perifocal frame.
"""
struct Perifocal <: AbstractFrame end

"""
Cartesian state which describes a spacecraft or body's position and velocity with respect to _something_.
"""
mutable struct CartesianState{F, LU, TU, FR} <: AbstractState{F, LU, TU, FR} 
    t::F
    r::SubArray{F, 1, MVector{6, F}, Tuple{UnitRange{Int64}}, true}
    v::SubArray{F, 1, MVector{6, F}, Tuple{UnitRange{Int64}}, true}
    rv::MVector{6,F}

    function CartesianState(r::AbstractVector{R}, v::AbstractVector{V}, epoch::E=0, frame=Bodycentric; lengthunit=u"km", timeunit=u"s") where {R<:Real, V<:Real, E<:Real}
        F  = promote_type(R, V, E)
        if !(F <: AbstractFloat)
            @warn "Type provided ($(string(F))) is not a float: defaulting to Float64."
            F = Float64
        end
        @assert length(r) == length(v) == 3 "Both arguments must have length equal to 3!"
        rv = MVector{6, F}(r[1], r[2], r[3], v[1], v[2], v[3])
        pos = @views rv[1:3]
        vel = @views rv[4:6]
        return new{F, typeof(lengthunit), typeof(timeunit), frame}(F(epoch), pos, vel, rv)
    end

    function CartesianState(r::AbstractVector{R}, v::AbstractVector{V}, epoch::Unitful.Time = 0u"s", frame=Bodycentric) where {R<:Unitful.Length, V<:Unitful.Velocity}

        # Get ready to commit a crime... we need the Time unit provided in velocity quantity V

        # V is a quantity with units Length / Time.  Let's reverse that so we have Time instead of Time^(-1)
        TL = inv(unit(V))

        # Now we need to find which index (1 or 2) contains the Time unit
        timeaxis = findfirst(T -> T isa Unitful.Dimension{:Time}, collect(typeof(typeof(TL).parameters[2]).parameters[1]))
        
        # Now that we have the proper index, let's select the time unit
        timeunit = Unitful.FreeUnits{(typeof(TL).parameters[1][timeaxis],), Unitful.𝐓, nothing}()

        # This is easy...
        lengthunit = unit(R)

        # Phew!
        return CartesianState(ustrip.(lengthunit, r), ustrip.(lengthunit/timeunit, v), ustrip(timeunit, epoch), frame; lengthunit = lengthunit, timeunit = timeunit)
    end

    function CartesianState(cart::CartesianState{F,LU,TU,FR}) where {F,LU,TU,FR}
        return CartesianState(cart.r,  cart.v, cart.t, FR; lengthunit = LU(), timeunit = TU())
    end

    function CartesianState(arr::StaticVector{6})
        return CartesianState(arr[1:3], arr[4:6])
    end
end

Base.getindex(cart::CartesianState, i::Int) = cart.rv[i]
Base.setindex!(cart::CartesianState, value, i::Int) = (cart.rv[i] = value)
Base.length(::CartesianState) = 6
Base.size(::CartesianState) = (6,)

function Base.convert(::Type{CartesianState{F,LU,TU}}, cart::CartesianState) where {F,LU,TU}
    r = F.(ustrip.(LU(), position_vector(cart)))
    v = F.(ustrip.(LU()/TU(), velocity_vector(cart)))
    t = F.(ustrip.(TU(), epoch(cart)))
    return CartesianState(r, v, t, coordinateframe(cart); lengthunit = LU(), timeunit = TU())
end

epoch(cart::CartesianState) = cart.t * timeunit(cart)


const NormalizedLengthUnit   = Unitful.FreeUnits{(), Unitful.𝐋, nothing}
const NormalizedTimeUnit     = Unitful.FreeUnits{(), Unitful.𝐓, nothing}

Base.eltype(::AbstractUnitfulStructure{F}) where F = F

"""
Returns the `Unitful.Length` unit associated with the state.
"""
lengthunit(::AbstractUnitfulStructure{F, LU, TU}) where {F,LU,TU} = LU()

"""
Returns the `Unitful.Time` unit associated with the state.
"""
timeunit(::AbstractUnitfulStructure{F, LU, TU}) where {F, LU, TU} = TU()

"""
Returns the `Unitful.Velocity` unit associated with the state.
"""
velocityunit(state::AbstractUnitfulStructure) = lengthunit(state) / timeunit(state)

"""
Returns the `MassParameter` unit associated with the state.
"""
massparameterunit(state::AbstractUnitfulStructure) = lengthunit(state)^3 / timeunit(state)^2

"""
Returns the coordinate frame.
"""
coordinateframe(::AbstractState{F, LU, TU, FR}) where {F, LU, TU, FR} = FR

function Base.show(io::IO, state::CartesianState{F,LU,TU,FR}) where {F,LU,TU,FR} 
    println(io, "  ", string(FR), " Cartesian State:")
    println("")
    println(io, "    t:  ", state.t, " ", string(TU()))
    println(io, "    r = ", [state.r[1] state.r[2] state.r[3]], " ", string(LU()))
    println(io, "    v = ", [state.v[1] state.v[2] state.v[3]], " ", string(LU()/TU()))
end

Base.show(io::IO, ::MIME"text/plain", state::CartesianState{F,LU,TU,FR}) where {F,LU,TU,FR} = show(io, state)

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