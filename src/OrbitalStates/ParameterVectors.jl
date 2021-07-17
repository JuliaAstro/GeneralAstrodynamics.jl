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
abstract type ParameterVector{F, MU, LU, TU, AU, N, T, B} <: ParameterizedLabelledArray{F, 1, T, SLArray{Tuple{N},F,1,N,T}} end

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
struct R2BPParameters{F, MU, LU, TU, AU, B} <: ParameterVector{F, MU, LU, TU, AU, 1, (:μ,), B}
    __rawdata::SLArray{Tuple{1}, F, 1, 1, (:μ,)}
    
    function R2BPParameters(μ; massunit=u"kg", lengthunit=u"km", timeunit=u"s", angularunit=u"°", name=:Primary)
        μμ = (μ isa MassParameter) ? ustrip(massunit^3/timeunit^2, μ) : μ
        F = eltype(μμ)
        F isa AbstractFloat || (F = Float64)
        B = Symbol(name)
        return new{F, massunit, lengthunit, timeunit, angularunit, B}(SLArray{Tuple{1}, F, 1, 1, (:μ,)}(μμ))
    end
end

"""
$(TYPEDEF)

All parameters required for the Circular Restricted Three-body Problem.
"""
struct CR3BPParameters{F, MU, LU, TU, AU, B} <: ParameterVector{F, MU, LU, TU, AU, 3, (:μ,:μ₁,:μ₂), B}
    __rawdata::SLArray{Tuple{3}, F, 1, 3, (:μ,:μ₁,:μ₂)}
    
    function CR3BPParameters(μ₁, μ₂; massunit=u"kg", lengthunit=u"km", timeunit=u"s", angularunit=u"°", primary=:Primary, secondary=:Secondary)
        μμ₁ = (μ₁ isa MassParameter) ? ustrip(massunit^3/timeunit^2, μ₁) : μ₁
        μμ₂ = (μ₂ isa MassParameter) ? ustrip(massunit^3/timeunit^2, μ₂) : μ₂
        μμ  = max(μμ₁, μμ₂) / (μμ₁+μμ₂)
        @assert μμ ≤ 1//2 "Nondimensional mass parameter must be less than 1/2, by definition! Received max(μ₁,μ₂)/(μ₁+μ₂)=$μμ."
        F = promote_type(typeof(μμ₁), typeof(μμ₂))
        F isa AbstractFloat || (F = Float64)
        B = (Symbol(primary), Symbol(secondary))
        return new{F, massunit, lengthunit, timeunit, angularunit, B}(SLArray{Tuple{3}, F, 1, 3, (:μ,:μ₁,:μ₂)}((max(μμ₁,μμ₂)/(μμ₁+μμ₂), μμ₁, μμ₂)))
    end

    function CR3BPParameters(μ::Number; massunit=u"kg", lengthunit=u"km", timeunit=u"s", angularunit=u"°", primary=:Unknown, secondary=:Unknown)
        μμ = (μ isa MassParameter) ? ustrip(massunit^3/timeunit^2, μ) : μ
        @assert μμ ≤ 1//2 "Nondimensional mass parameter must be less than 1/2, by definition! Received $μμ."
        F = typeof(μμ)
        F isa AbstractFloat || (F = Float64)
        B = (Symbol(primary), Symbol(secondary))
        return new{F, massunit, lengthunit, timeunit, angularunit, B}(SLArray{Tuple{3}, F, 1, 3, (:μ,:μ₁,:μ₂)}(μμ, NaN , NaN))
    end
end

