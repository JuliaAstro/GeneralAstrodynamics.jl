using Core: Argument
#
# Provides `create_array` methods for `StateVector` and `ParameterVector`
# instances!
#

@inline function SymbolicUtils.Code.create_array(A::Type{<:StateVector}, T, nd::Val, d::Val{dims}, elems...) where {dims}
    data = SymbolicUtils.Code.create_array(SArray, T, nd, d, elems...)
    return (A)(data.data)
end
