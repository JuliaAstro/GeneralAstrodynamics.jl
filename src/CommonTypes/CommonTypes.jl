"""
Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module CommonTypes

using Reexport
@reexport using Unitful, UnitfulAngles, UnitfulAstro

include("../Misc/DocStringExtensions.jl")
include("../Misc/UnitfulAliases.jl")

export AbstractBody, AbstractOrbitalSystem, AbstractTrajectory, AbstractState

""" 
Abstract type for bodies in space: both `CelestialBody`s (in
`TwoBody.jl`), and `Body`s (in `NBody.jl`).
"""
abstract type AbstractBody end

"""
Abstract type describing select Astrodynamics problems.
"""
abstract type AbstractOrbitalSystem end

"""
Abstract type describing a collection of states resulting from numerical integration
"""
abstract type AbstractTrajectory end

using StaticArrays, Unitful, ComponentArrays

"""
Abstract type describing an orbital state.
"""
abstract type AbstractState end

"""
Cartesian state which describes a spacecraft or body's position and velocity with respect to _something_.
"""
struct CartesianState{F<:AbstractFloat, LU, TU} <: AbstractState where {LU<:Unitful.LengthUnits, TU<:Unitful.TimeUnits}
    state::ComponentArray{F,1,Vector{F},Tuple{Axis{(r = 1:3, v = 4:6)}}}
    
    function CartesianState(r::AbstractVector{R}, v::AbstractVector{V}; lengthunit=u"km", timeunit=u"s") where {R<:Real, V<:Real}
        F  = promote_type(R, V)
        if !(F <: AbstractFloat)
            @warn "Type provided ($(string(F))) is not a float: defaulting to Float64."
            F = Float64
        end
        return new{F, lengthunit, timeunit}(ComponentVector{F}(r=r, v=v))
    end

    function CartesianState(state::ComponentArray{F,1,Vector{F},Tuple{Axis{(r = 1:3, v = 4:6)}}}; lengthunit=u"km", timeunit=u"s") where F<:Real
        if !(F <: AbstractFloat)
            @warn "Type provided ($(string(F))) is not a float: defaulting to Float64."
            return new{Float64, lengthunit, timeunit}(Float64.(state))
        else
            return new{F, lengthunit, timeunit}(F.(state))
        end
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

CartesianState(ones(3) * u"km", zeros(3) * u"km/s")