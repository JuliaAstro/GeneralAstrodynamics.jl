#
# Provides concrete types under `ParameterizedLabelledArray`
# which contain parameters to describe astrodynamic systems.
#

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
massparamunit(::ParameterVector{F, MU, LU, TU, AU}) where {F, MU, LU, TU, AU} = LU^3/TU^2

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
name(::ParameterVector{F, MU, LU, TU, AU, N, T, B}) where {F, MU, LU, TU, AU, N, T, B} = B

"""
All parameters required for the Restricted Two-body Problem.
"""
struct R2BPParameters{F, MU, LU, TU, AU, B} <: ParameterVector{F, MU, LU, TU, AU, 1, (:μ,), B}
    __rawdata::SLArray{Tuple{1}, F, 1, 1, (:μ,)}
    
    function R2BPParameters(μ; massunit=u"kg", lengthunit=u"km", timeunit=u"s", angularunit=u"°", name=:Primary)
        μμ = (μ isa Unitful.AbstractQuantity) ? ustrip(lengthunit^3/timeunit^2, μ) : μ
        F = eltype(μμ)
        F isa AbstractFloat || (F = Float64)
        B = Symbol(name)
        return new{F, massunit, lengthunit, timeunit, angularunit, B}(SLArray{Tuple{1}, F, 1, 1, (:μ,)}(μμ))
    end
end

"""
Returns the normalized mass parameter, `μ`.
"""
get_μ(sys::R2BPParameters) = sys[1]

"""
Displays `R2BPParameters`.
"""
function Base.show(io::IO, state::R2BPParameters; showfloats=true, space="")
    println(io, space, "R2BP parameters ", name(state) == "Unknown" ? "" : "for the $(name(state)) system", showfloats ? " with eltype $(eltype(state))" : "", "\n")
    println(io, space, "  ", "μ = ", get_μ(state), " ", massparamunit(state))
end
Base.show(io::IO, ::MIME"text/plain", state::R2BPParameters; kwargs...) = Base.show(io, state; kwargs...)

"""
Converts types and units for a `R2BPParameters`.
"""
function Base.convert(::Type{R2BPParameters{F, MU, LU, TU, AU}}, system::R2BPParameters) where {F, MU, LU, TU, AU}
    μ = ustrip(LU^3/TU^2, get_μ(system) * massparamunit(system)) |> F
    return R2BPParameters(μ; massunit=MU, lengthunit=LU, timeunit=TU, angularunit=AU) 
end

"""
All parameters required for the Circular Restricted Three-body Problem.
"""
struct CR3BPParameters{F, MU, LU, TU, AU, B} <: ParameterVector{F, MU, LU, TU, AU, 1, (:μ,), B}
    __rawdata::SLArray{Tuple{1}, F, 1, 1, (:μ,)}
    
    function CR3BPParameters(μ::Real; massunit=u"kg", lengthunit=missing, timeunit=missing, angularunit=u"°", primary=:Primary, secondary=:Secondary)
        @assert μ ≤ 1//2 "Nondimensional mass parameter must be less than 1/2, by definition! Received $μ."
        F = typeof(μ)
        F isa AbstractFloat || (F = Float64)
        B = (Symbol(primary), Symbol(secondary))
        return new{F, massunit, lengthunit, timeunit, angularunit, B}(SLArray{Tuple{1}, F, 1, 1, (:μ,)}(μ))
    end

    function CR3BPParameters(μ₁::Unitful.AbstractQuantity, μ₂::Unitful.AbstractQuantity, a::Unitful.Length; massunit=u"kg", angularunit=u"°", primary=:Primary, secondary=:Secondary)
        ∑μᵢ = μ₁ + μ₂
        μ  = min(μ₁, μ₂) / ∑μᵢ
        TU = 2π * upreferred(√(a^3 / ∑μᵢ))
        F = typeof(μ)
        F isa AbstractFloat || (F = Float64)
        B = (Symbol(primary), Symbol(secondary))
        return new{F, massunit, a, TU, angularunit, B}(SLArray{Tuple{1}, F, 1, 1, (:μ)}(F(μ)))
    end
end

"""
Returns the normalized mass parameter, `μ`.
"""
get_μ(sys::CR3BPParameters) = sys[1]

"""
Displays a `CR3BPParameters` instance.
"""
function Base.show(io::IO, sys::CR3BPParameters; showfloats=true, space="") 
    sysname   = all(name(sys) .== (:Primary, :Secondary)) ? "" : begin (n1, n2) = name(sys); string(" for the ", n1, "-", n2, " system") end
    println(io, space, "CR3BP parameters", sysname, showfloats ? " with eltype $(eltype(sys))" : "", "\n")
    println(io, space, "  ", "μ = ", get_μ(sys))
end
Base.show(io::IO, ::MIME"text/plain", sys::CR3BPParameters; kwargs...) = Base.show(io::IO, sys::CR3BPParameters; kwargs...)

"""
Converts types and units for a `CR3BPParameters`.
"""
function Base.convert(::Type{CR3BPParameters{F, MU, LU, TU, AU}}, system::CR3BPParameters) where {F, MU, LU, TU, AU}
    μ  = get_μ(system)
    μ₁, μ₂ = let 
        ∑μᵢ = TU^3 / (TU / 2π)^2
        μ₂ = μ * ∑μᵢ
        μ₁ = ∑μᵢ - μ₂

        μ₁, μ₂
    end

    μₙ = min(μ₁, μ₂) / (μ₁ + μ₂)
    return CR3BPParameters(μₙ; massunit=MU, lengthunit=LU, timeunit=TU, angularunit=AU)
end