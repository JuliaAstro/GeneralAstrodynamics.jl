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
A supertype for parameter representations in astrodynamics.
"""
abstract type ParameterVector{F, MU, LU, TU, AU, N, T, B} <: ParameterizedLabelledArray{F, 1, T, SLArray{Tuple{N},F,1,N,T}} end

"""
Returns the mass unit of the parameter vector.
"""
Base.@pure massunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = MU


"""
Returns the length unit of the parameter vector.
"""
Base.@pure lengthunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = LU

"""
Returns the time unit of the parameter vector.
"""
Base.@pure timeunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = TU

"""
Returns the angular unit of the parameter vector.
"""
Base.@pure angularunit(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = AU

"""
Returns the velocity unit of the state vector.
"""
velocityunit(::ParameterVector{F, MU, LU, TU, AU}) where {F, MU, LU, TU, AU} = LU/TU

"""
Returns the mass-parameter unit of the state vector.
"""
massparamunit(::ParameterVector{F, MU, LU, TU, AU}) where {F, MU, LU, TU, AU} = MU^3/LU^2

"""
Returns the length of the parameter vector.
"""
Base.@pure Base.length(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = N

"""
Returns the size of the parameter vector.
"""
Base.@pure Base.size(::ParameterVector{F, MU, LU, TU, AU, N}) where {F, MU, LU, TU, AU, N} = (N,)


"""
Returns the name of the parameter vector.
"""
name(::ParameterVector{F, MU, LU, TU, AU, N, T, B}) where {F, MU, LU, TU, AU, N, T, B} = string(B)

"""
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
Converts between `R2BPParameters`.
"""
Base.convert(::Type{R2BPParameters{F, MU, LU, TU, AU}}, params::R2BPParameters) where {F,MU,LU,TU,AU} = R2BPParameters(
    F(get_μ(params)); massparamunit = massparamunit(params), 
                      lengthunit    = lengthunit(params),
                      timeunit      = timeunit(params),
                      angularunit   = angularunit(params),
                      name          = name(params)
)

"""
Returns the mass parameter of the R2BP system.
"""
massparameter(system::R2BPParameters) = system[1] * massparamunit(system)

"""
Displays `R2BPParameters`.
"""
function Base.show(io::IO, state::R2BPParameters; showfloats=true, space="")
    println(io, space, "R2BP Parameters ", name(state) == "Unknown" ? "" : "for the $(name(state)) System", showfloats ? " ($(eltype(state)))" : "")
    println(io, space, "  μ = $(state[1]) $(lengthunit(state)^3/timeunit(state)^2)")
end

"""
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
        return new{F, massunit, lengthunit, timeunit, angularunit, B}(SLArray{Tuple{3}, F, 1, 3, (:μ,:μ₁,:μ₂)}((max(μμ₁,μμ₂)/(μμ₁+μμ₂), max(μμ₁, μμ₂), min(μμ₁, μμ₂))))
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

"""
Converts between `CR3BPParameters`.
"""
function Base.convert(::Type{CR3BPParameters{F, MU, LU, TU, AU}}, params::CR3BPParameters) where {F, MU, LU, TU, AU}
    if any(isnan, (params.μ₁, params.μ₂))
        return CR3BPParameters(
            F(get_μ(params)); massparamunit = massparamunit(params), 
            lengthunit    = lengthunit(params),
            timeunit      = timeunit(params),
            angularunit   = angularunit(params),
            name          = name(params)
        )
    else 
        return CR3BPParameters(
            F(get_μ₁(params)), F(get_μ₂(params)); 
            massparamunit = massparamunit(params), 
            lengthunit    = lengthunit(params),
            timeunit      = timeunit(params),
            angularunit   = angularunit(params),
            name          = name(params)
        )
    end
end

"""
Displays `CR3BPParameters`.
"""
function Base.show(io::IO, state::CR3BPParameters; showfloats=true, space="")
    println(io, space, "CR3BP Parameters ", any(x == "Unknown", name(state)) ? "" : "for the $(name(state)[1])-$(name(state)[2]) System", showfloats ? " ($(eltype(state)))" : "")
    println(io, space, "  μ = $(state[1]) $(lengthunit(state)^3/timeunit(state)^2)")
end


"""
Returns the mass parameter of the CR3BP system.
"""
normalized_massparameter(system::CR3BPParameters) = system[1]

"""
Returns the mass parameters of the CR3BP system.
"""
massparameters(system::CR3BPParameters) = (system[2], system[3]) * massparamunit(system)

"""
Returns the primary mass parameter of the CR3BP system.
"""
primary_massparameter(system::CR3BPParameters) = system[2] * massparamunit(system)

"""
Returns the secondary mass parameter of the CR3BP system.
"""
secondary_massparameter(system::CR3BPParameters) = system[3] * massparamunit(system)