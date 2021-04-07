"""
Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module CommonTypes


macro boilerplate(struct_definition)
    firstline_post_struct = split(split(string(struct_definition), "\n")[1], "struct ")[2]
    if ' ' ∈ firstline_post_struct
        structname = split(firstline_post_struct, " ")  |> first
    end
    if '{' ∈ structname
        structname = split(structname, "{") |> first
    end

    if '\n' ∈ structname
        structname = split(firstline_post_struct, "\n") |> first
    end

    structname = structname |> Symbol

    convert = :(Base.convert(::Type{T}, o::$(esc(structname))) where {T<:Number} = $(esc(structname))(map(x-> typeof(x) <: Union{AbstractArray{<:Number}, Number} ? T.(x) : x, fieldnames($(esc(structname))))...))
    if '{' ∈ firstline_post_struct
        promote = :(Base.promote(::Type{$(esc(structname)){A}}, ::Type{$(esc(structname)){B}}) where {A,B} = $(esc(structname)){promote_type(A,B)})
    else
        promote = :(Base.promote(::Type($(esc(structname)), ::Type{$(esc(structname))})) = $(esc(structname)))
    end
    Float16  = :(Core.Float16(o::$(esc(structname))) = convert(Float16, o))
    Float32  = :(Core.Float32(o::$(esc(structname))) = convert(Float32, o))
    Float64  = :(Core.Float64(o::$(esc(structname))) = convert(Float64, o))
    BigFloat = :(Base.MPFR.BigFloat(o::$(esc(structname))) = convert(BigFloat, o))
        
    quote
        $struct_definition
        $convert
        $promote
        $Float16
        $Float32
        $Float64
    end 
end

macro export_boilerplate()
    return :(export convert, promote, Float16, Float32, Float64)
end

export AbstractBody, AbstractOrbitalSystem, AbstractTrajectory, AbstractCartesianState, CartesianState
export getindex, setindex!, lengthunit, timeunit, velocityunit
export position_vector, velocity_vector, scalar_position, scalar_velocity
# export @boilerplate, @export_boilerplate
# @export_boilerplate

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
abstract type AbstractOrbitalSystem end

"""
Abstract type describing a collection of states resulting from numerical integration
"""
abstract type AbstractTrajectory end

"""
Abstract type describing an orbital state.
"""
abstract type AbstractCartesianState{F<:AbstractFloat} <: StaticVector{6,F} end

mutable struct CartesianState{F<:AbstractFloat, LU, TU} <: AbstractCartesianState{F} where {LU<:Unitful.LengthFreeUnits, TU<:Unitful.TimeFreeUnits}
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

@doc "Cartesian state which describes a spacecraft or body's position and velocity with respect to _something_." CartesianState

# We need beyond the boilerplate for `CartesianState`...
Base.convert(::Type{T}, o::CartesianState) where {T<:CartesianState} = CartesianState(o.r / upreferred((1 * T.parameters[2]) / (1 * lengthunit(o))), o.v / upreferred((1 * T.parameters[2]) / (1 * lengthunit(o))) * upreferred((1 * T.parameters[3]) * (1 * timeunit(o))); lengthunit = T.parameters[2], timeunit = T.parameters[3])

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

function Base.show(io::IO, ::MIME"text/plain", state::CartesianState{F, LU, TU}) where {F, LU, TU} 
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