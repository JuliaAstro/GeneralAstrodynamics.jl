#
# Provides concrete types under `ParameterizedLabelledArray`
# which contain parameters to describe astrodynamic systems.
#

    
Unitful.@derived_dimension MassParameter Unitful.𝐋^3/Unitful.𝐓^2

@doc """
A unit dimension alias for length^3 / time^2. This is a common
dimension used in astrodynamics calculations.
"""
MassParameter

"""
$(TYPEDEF)

A supertype for parameter representations in astrodynamics.
"""
abstract type ParameterVector{F, MU, LU, TU, AU, N, T} <: ParameterizedLabelledArray{F, 1, T, SLArray{Tuple{N},F,1,N,T}} end

"""
$(SIGNATURES)

Returns the mass unit of the parameter vector.
"""
Base.@pure massunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = MU


"""
$(SIGNATURES)

Returns the length unit of the parameter vector.
"""
Base.@pure lengthunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = LU

"""
$(SIGNATURES)

Returns the time unit of the parameter vector.
"""
Base.@pure timeunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = TU

"""
$(SIGNATURES)

Returns the angular unit of the parameter vector.
"""
Base.@pure angularunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = AU

"""
$(SIGNATURES)

Returns the length of the parameter vector.
"""
Base.@pure Base.length(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = N

"""
$(SIGNATURES)

Returns the size of the parameter vector.
"""
Base.@pure Base.size(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = (N,)

"""
$(TYPEDEF)

All parameters required for the Restricted Two-body Problem.
"""
struct R2BPParameters{F, MU, LU, TU, AU} <: ParameterVector{F, MU, LU, TU, AU, 1, (:μ,)}
    __rawdata::SLArray{Tuple{1}, F, 1, 1, (:μ,)}
    
    function R2BPParameters(μ; massunit=u"kg", lengthunit=u"km", timeunit=u"s", angularunit=u"°")
        μμ = (μ isa MassParameter) ? ustrip(massunit^3/timeunit^2, μ) : μ
        F = eltype(μμ)
        F isa AbstractFloat || (F = Float64)
        return new{F, massunit, lengthunit, timeunit, angularunit}(SLArray{Tuple{1}, F, 1, 1, (:μ,)}(μμ))
    end
end

