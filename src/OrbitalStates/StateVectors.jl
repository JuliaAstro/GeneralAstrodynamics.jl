#
# State vector descriptions.
#

"""
$(TYPEDEF)

A supertype for all state representations in astrodynamics.
"""
abstract type StateVector{F, LU, TU, AU} <: DenseVector{F} end

"""
$(SIGNATURES)

Returns the lengthunit of the state vector.
"""
Base.@pure lengthunit(::StateVector{F, LU, TU, AU}) where {F, LU, TU, AU} = LU

"""
$(SIGNATURES)

Returns the timeunit of the state vector.
"""
Base.@pure timeunit(::StateVector{F, LU, TU, AU}) where {F, LU, TU, AU} = TU

"""
$(SIGNATURES)

Returns the angularunit of the state vector.
"""
Base.@pure angularunit(::StateVector{F, LU, TU, AU}) where {F, LU, TU, AU} = AU

"""
$(SIGNATURES)

The length of any `StateVector` is 6!
"""
Base.@pure Base.length(::StateVector) = 6

"""
$(SIGNATURES)

The size of any `StateVector` is (6,)!
"""
Base.@pure Base.size(::StateVector) = (6,)

"""
$(SIGNATURES)

Returns the index of the `StateVector`.
"""
Base.getindex(state::StateVector, args...) = Base.getindex(state.__rawdata, args...)

"""
$(SIGNATURES)

Sets the index of the `StateVector`.
"""
Base.setindex!(state::StateVector, args...) = Base.setindex!(state.__rawdata, args...)

"""
$(SIGNATURES)

Returns the memory stride for any `StateVector`.
"""
Base.elsize(::StateVector{F}) where F = sizeof(F)

#=
All code in the following begin ... end block was copied, and modified 
from LabelledArrays.jl source code. Their license is shown below. 
All credit goes to LabelledArrays.jl developers.

The LabelledArrays.jl package is licensed under the MIT "Expat" License:

> Copyright (c) 2017: Christopher Rackauckas.
>
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
>
> of this software and associated documentation files (the "Software"), to deal
>
> in the Software without restriction, including without limitation the rights
>
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
>
> copies of the Software, and to permit persons to whom the Software is
>
> furnished to do so, subject to the following conditions:
>
>
>
> The above copyright notice and this permission notice shall be included in all
>
> copies or substantial portions of the Software.
>
>
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
>
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
>
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
>
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
>
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
>
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
>
> SOFTWARE.
>
>
=#
begin
    """
    $(SIGNATURES)

    Overrides `Base.getproperty` for all `StateVector` instances.
    """
    Base.@propagate_inbounds function Base.getproperty(x::StateVector, s::Symbol)
        if s ∈ propertynames(x)
            return getproperty(getfield(x, :__rawdata), s)
        end
        return getfield(x, s) # will throw an error if s is not :__rawdata!
    end

    """
    $(SIGNATURES)

    Sets indices of a `StateVector` via label.
    """
    Base.@propagate_inbounds function Base.setproperty!(x::StateVector, s::Symbol,y)
        if s ∈ propertynames(x)
            return setproperty!(getfield(x, :__rawdata), s, y)
        end
        setfield!(x, s, y)
    end

    """
    $(SIGNATURES)

    Overrides `similar` for `StateVector` instances.
    """
    function Base.similar(x::StateVector, ::Type{S}, dims::NTuple{N,Int}) where {S,N}
        tmp = similar(x.__rawdata, S, dims)
        return typeof(x)(tmp)
    end
    
    """
    $(SIGNATURES)

    Shallow copies a `StateVector`.
    """
    Base.copy(x::StateVector) = typeof(x)(copy(getfield(x,:__rawdata)))

    """
    $(SIGNATURES)

    Deep copies a `StateVector`.
    """
    Base.deepcopy(x::StateVector) = typeof(x)(deepcopy(getfield(x, :__rawdata)))

    """
    $(SIGNATURES)

    Copies one `StateVector` to another.
    """
    Base.copyto!(x::C,y::C) where C <: StateVector = copyto!(getfield(x,:__rawdata),getfield(y,:__rawdata))
    
    """
    $(SIGNATURES)

    Provides `unsafe_convert` for `StateVector` types for use with LAPACK.
    """
    Base.unsafe_convert(::Type{Ptr{T}}, a::StateVector{F}) where {T, F} = Base.unsafe_convert(Ptr{T}, getfield(a,:__rawdata))
    
    """
    $(SIGNATURES)

    Converts the underlying floating point type for a `StateVector`.
    """
    Base.convert(::Type{T}, x) where {T<:StateVector} = T(x)

    """
    $(SIGNATURES)

    Converts the underlying floating point type for a `StateVector`.
    """
    Base.convert(::Type{T}, x::T) where {T<:StateVector} = x

    """
    $(SIGNATURES)

    Converts the underlying floating point type for a `StateVector`.
    """
    Base.convert(::Type{<:Array},x::StateVector) = convert(Array, getfield(x,:__rawdata))
    
    """
    $(SIGNATURES)

    Reshapes a `StateVector`.
    """
    ArrayInterface.restructure(x::StateVector{F1}, y::StateVector{F2}) where {F1, F2} = reshape(y, size(x)...)
    
    """
    $(SIGNATURES)

    Implements `dataids` for a `StateVector` instance.
    """
    Base.dataids(A::StateVector) = Base.dataids(A.__rawdata)
end

"""
$(TYPEDEF)

A Cartesian state vector with length 6. Internally
uses `MVector` and `LVector` to store data.
Data is accessible via labels, which are
`x, y, z, ẋ, ẏ, ż, r, v`, where `r` 
access `x,y,z` and `v` accesses `ẋ, ẏ, ż`.
"""
mutable struct CartesianState{F, LU, TU, AU} <: StateVector{F, LU, TU, AU} 
    __rawdata::L where L<:LArray{F, 1, MVector{6, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)}

    function CartesianState(statevector; lengthunit=u"km", timeunit=u"s", angularunit=u"°")
        F = eltype(statevector)
        F isa AbstractFloat || (F = Float64)

        return new{F, lengthunit, timeunit, angularunit}(
            statevector
        )
    end
end

"""
$(SIGNATURES)

Constructs a `CartesianState` from any `AbstractVector`.
"""
CartesianState{F, LU, TU, AU}(statevector::AbstractVector{<:Real}) where {F, LU, TU, AU} = CartesianState(
    convert(LArray{F, 1, MVector{6, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)}, statevector)
)

"""
$(SIGNATURES)

Constructs a `CartesianState` from provided position and velocity vectors.
"""
function CartesianState(r::AbstractVector, v::AbstractVector; 
                        lengthunit=(unit(eltype(r)) isa Unitful.Length ? unit(eltype(r)) : u"km"),
                        timeunit=(unit(eltype(v)) isa Unitful.Velocity ? lengthunit / unit(eltype(v)) : u"s"),
                        angularunit=u"°")
    rr = (eltype(r) isa Unitful.Length ? ustrip.(lengthunit, r) : r)
    vv = (eltype(v) isa Unitful.Velocity ? ustrip.(lengthunit / timeunit, v) : v)
    F = promote_type(eltype(rr), eltype(vv))
    F isa AbstractFloat || (F = Float64)
    return CartesianState{F, lengthunit, timeunit, angularunit}(vcat(rr,vv))
end

"""
$(SIGNATURES)

Returns dot-accessible property names for a `CartesianState`: (:x, :y, :z, :ẋ, :ẏ, :ż, :r, :v).
"""
Base.@pure Base.propertynames(::CartesianState) =  (:x, :y, :z, :ẋ, :ẏ, :ż, :r, :v)

"""
$(TYPEDEF)

A Keplerian state vector with length 6. Internally
uses `MVector` and `LVector` to store data.
Data is accessible via labels, which are
`e, a, i, Ω, ω, ν`.
"""
mutable struct KeplerianState{F, LU, TU, AU} <: StateVector{F, LU, TU, AU} 
    __rawdata::L where L<:LArray{F, 1, MVector{6, F}, (e=1, a=2, i=3, Ω=4, ω=5, ν=6)}

    function KeplerianState(statevector; lengthunit=u"km", timeunit=u"s", angularunit=u"°")
        F = eltype(statevector)
        F isa AbstractFloat || (F = Float64)

        return new{F, lengthunit, timeunit, angularunit}(
            statevector
        )
    end
end

"""
$(SIGNATURES)

Constructs a `KeplerianState` from any `AbstractVector`.
"""
KeplerianState{F, LU, TU, AU}(statevector::AbstractVector{<:Real}) where {F, LU, TU, AU} = KeplerianState(
    convert(LArray{F, 1, MVector{6, F}, (e=1, a=2, i=3, Ω=4, ω=5, ν=6)}, statevector)
)

"""
$(SIGNATURES)

Constructs a `KeplerianState` from provided position and velocity vectors.
"""
function KeplerianState(e::Real, a, i, Ω, ω, ν; 
                        lengthunit=(a isa Unitful.Length ? unit(a) : u"km"),
                        timeunit=u"s",
                        angularunit=(i isa DimensionlessQuantity ? unit(i) : u"°"))
    aa = (a isa Unitful.Length ? ustrip(lengthunit, a) : a)
    ii = (i isa Unitful.DimensionlessQuantity ? ustrip(angularunit, i) : i)
    ΩΩ = (Ω isa Unitful.DimensionlessQuantity ? ustrip(angularunit, Ω) : Ω)
    ωω = (ω isa Unitful.DimensionlessQuantity ? ustrip(angularunit, ω) : ω)
    νν = (ν isa Unitful.DimensionlessQuantity ? ustrip(angularunit, ν) : ν)
    F = promote_type(typeof(e), typeof(aa), typeof(ii), typeof(ΩΩ), typeof(ωω), typeof(νν))
    F isa AbstractFloat || (F = Float64)
    return Keplerian{F, lengthunit, timeunit, angularunit}(@LArray(MVector{6,F}(e,a,i,Ω,ω,ν),(e=1, a=2, i=3, Ω=4, ω=5, ν=6)))
end

"""
$(SIGNATURES)

Returns dot-accessible property names for a `CartesianState`: (:x, :y, :z, :ẋ, :ẏ, :ż, :r, :v).
"""
Base.@pure Base.propertynames(::KeplerianState) = (:e, :a, :i, :Ω, :ω, :ν)
