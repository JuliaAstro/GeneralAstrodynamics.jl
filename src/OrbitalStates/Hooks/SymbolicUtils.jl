#
# Provides `create_array` methods for `StateVector` and `ParameterVector`
# instances!
#

@inline function create_array(A::Type{<:StateVector}, a::Nothing, d::Val{dims}, elems...) where dims
    return A((elems...,))
end

# and

@inline function create_array(::Type{<:MyArray}, T, ::Val{dims}, elems...) where dims
    return A(T.(elems...,))
end