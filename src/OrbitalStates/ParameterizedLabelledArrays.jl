#
# Using LabelledArrays source code and types to 
# "pass through" all properties to an underlying
# LArray or SLArray. This allows us to write 
# a new abstract type that *functions*
# like a LabelledArray type!
#
# All code in this file was copied, and modified 
# from LabelledArrays.jl source code. Their license is shown below. 
# All credit goes to LabelledArrays.jl developers.
#

const __LABELLED_ARRAYS_LICENSE = """
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
"""

const __LABELLED_ARRAYS_CREDITS = """
This source code which provides this functionality was copied directly from LabelledArrays.jl source code. The LabelledArrays.jl license text is provided in Julia's **Extended Help** section (accessible via `@doc`, or `??` in Julia's REPL).

# Extended Help

__LabelledArrays.jl License__

$__LABELLED_ARRAYS_LICENSE
"""

"""
$(TYPEDEF)

A supertype for types that *function* like 
`LabelledArray.LArray` or `LabelledArray.SLArray` 
instances, but are under a new type tree. This is
used in `GeneralAstrodynamics` for parameterizing 
astrodynamics state vectors by physical units.

!!! note
    All subtypes __must__ have only one field:
    a `LabelledArrays.LArray` or `LabelledArrays.SLArray`
    field called `__rawdata`. All methods on this abstract
    type require this field to be called `__rawdata`!

Nearly all code which acts on this type is copied 
and / or modified from LabelledArrays.jl source code.
All credit goes to LabelledArrays.jl developers. 
The LabelledArrays.jl LICENSE file is provided 
in this docstring under Julia's __Extended Help__
docstring section.

# Extended help

__LabelledArrays.jl License__

$__LABELLED_ARRAYS_LICENSE
"""
abstract type ParameterizedLabelledArray{F,N,T,L} <: DenseArray{F,N} end

"""
$(SIGNATURES)

Returns dot-accessible property names for a `ParameterizedLabelledArray`.
"""
Base.propertynames(::ParameterizedLabelledArray{F,N,T}) where {F,N,T} = keys(T)


"""
$(SIGNATURES)

Overrides `Base.getproperty` for all `ParameterizedLabelledArray` instances.

$__LABELLED_ARRAYS_CREDITS
"""
Base.@propagate_inbounds function Base.getproperty(x::ParameterizedLabelledArray, s::Symbol)
    if s ∈ propertynames(x)
        return getproperty(getfield(x, :__rawdata), s)
    end
    return getfield(x, s) # will throw an error if s is not :__rawdata!
end

"""
$(SIGNATURES)

Sets indices of a `ParameterizedLabelledArray` via label.

$__LABELLED_ARRAYS_CREDITS
"""
Base.@propagate_inbounds function Base.setproperty!(x::ParameterizedLabelledArray, s::Symbol,y)
    if s ∈ propertynames(x)
        return setproperty!(getfield(x, :__rawdata), s, y)
    end
    setfield!(x, s, y)
end

"""
$(SIGNATURES)

Overrides `similar` for `ParameterizedLabelledArray` instances.

$__LABELLED_ARRAYS_CREDITS
"""
function Base.similar(x::ParameterizedLabelledArray, ::Type{S}, dims::NTuple{N,Int}) where {S,N}
    tmp = similar(x.__rawdata, S, dims)
    return typeof(x)(tmp)
end

"""
$(SIGNATURES)

Shallow copies a `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.copy(x::ParameterizedLabelledArray) = typeof(x)(copy(getfield(x,:__rawdata)))

"""
$(SIGNATURES)

Deep copies a `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.deepcopy(x::ParameterizedLabelledArray) = typeof(x)(deepcopy(getfield(x, :__rawdata)))

"""
$(SIGNATURES)

Copies one `ParameterizedLabelledArray` to another.

$__LABELLED_ARRAYS_CREDITS
"""
Base.copyto!(x::C,y::C) where C <: ParameterizedLabelledArray = copyto!(getfield(x,:__rawdata),getfield(y,:__rawdata))

"""
$(SIGNATURES)

Provides `unsafe_convert` for `ParameterizedLabelledArray` types for use with LAPACK.

$__LABELLED_ARRAYS_CREDITS
"""
Base.unsafe_convert(::Type{Ptr{T}}, a::ParameterizedLabelledArray{F}) where {T, F} = Base.unsafe_convert(Ptr{T}, getfield(a,:__rawdata))

"""
$(SIGNATURES)

Converts the underlying floating point type for a `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.convert(::Type{T}, x) where {T<:ParameterizedLabelledArray} = T(x)

"""
$(SIGNATURES)

Converts the underlying floating point type for a `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.convert(::Type{T}, x::T) where {T<:ParameterizedLabelledArray} = x

"""
$(SIGNATURES)

Converts the underlying floating point type for a `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.convert(::Type{<:Array},x::ParameterizedLabelledArray) = convert(Array, getfield(x,:__rawdata))

"""
$(SIGNATURES)

Reshapes a `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
ArrayInterface.restructure(x::ParameterizedLabelledArray{F1}, y::ParameterizedLabelledArray{F2}) where {F1, F2} = reshape(y, size(x)...)

"""
$(SIGNATURES)

Implements `dataids` for a `ParameterizedLabelledArray` instance.

$__LABELLED_ARRAYS_CREDITS
"""
Base.dataids(A::ParameterizedLabelledArray) = Base.dataids(A.__rawdata)

"""
$(SIGNATURES)

Returns the index of the `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.getindex(state::ParameterizedLabelledArray, args...) = Base.getindex(state.__rawdata, args...)

"""
$(SIGNATURES)

Sets the index of the `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.setindex!(state::ParameterizedLabelledArray, args...) = Base.setindex!(state.__rawdata, args...)

"""
$(SIGNATURES)

Returns the memory stride for any `ParameterizedLabelledArray`.

$__LABELLED_ARRAYS_CREDITS
"""
Base.elsize(::ParameterizedLabelledArray{F}) where F = sizeof(F)