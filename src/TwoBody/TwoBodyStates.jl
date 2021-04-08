#
#   RestrictedTwoBodyStates.jl
#
#   Describes Two Body Orbits through Cartesian coordinates and Orbital Elements.
# 

"""
Abstract type for all four conic sections.
"""
abstract type AbstractConic end

"""
Type for orbits in the circular conic section.
"""
struct Circular <: AbstractConic end

"""
Type for orbits in the elliptical conic section.
"""
struct Elliptical <: AbstractConic end

"""
Type for orbits in the parabolic conic section.
"""
struct Parabolic <: AbstractConic end

"""
Type for orbits in the hyperbolic conic section.
"""
struct Hyperbolic <: AbstractConic end

"""
Type for invalid orbits (orbits with NaN fields)
"""
struct Invalid <: AbstractConic end

"""
Type representing large bodies in space. Currently, the following
solar system bodies are supported:

Sun, Mercury, Venus, Earth, Moon (Luna), Mars, Jupiter, 
Saturn, Uranus, Neptune, Pluto.
"""
struct RestrictedTwoBodySystem{F, LU, TU} <: AbstractUnitfulState{F, LU, TU}
    μ::F
    name::String

    function RestrictedTwoBodySystem(μ::Real, name=""; lengthunit = u"km", timeunit = u"s")
        T = typeof(μ)
        if !(T <: AbstractFloat)
            @warn "Non-float type $(string(T)) provided. Defaulting to Float64."
            T = Float64
        end
        return new{T, lengthunit, timeunit}(T(μ), name)
    end

    function RestrictedTwoBodySystem(μ::MassParameter, name="") 
        T = typeof(ustrip(μ))
        if !(T <: AbstractFloat)
            @warn "Non-float type $(string(T)) provided. Defaulting to Float64."
            T = Float64
        end

        M = unit(μ)
        lengthaxis = findfirst(T -> T isa Unitful.Dimension{:Length}, collect(typeof(typeof(M).parameters[2]).parameters[1]))
        timeaxis   = findfirst(T -> T isa Unitful.Dimension{:Time}, collect(typeof(typeof(M).parameters[2]).parameters[1]))

        lengthunit = (Unitful.FreeUnits{(typeof(M).parameters[1][lengthaxis], ), Unitful.𝐋^3}())^(1//3)
        timeunit   = (Unitful.FreeUnits{(typeof(M).parameters[1][timeaxis],), Unitful.𝐓^2}())^(-1//2)

        return new{T, lengthunit, timeunit}(T(ustrip(lengthunit^3 / timeunit^2, μ)), name)
    end
end

"""
Custom display for `RestrictedTwoBodySystem` instances.
"""
function Base.show(io::IO, body::RestrictedTwoBodySystem{F,LU,TU}) where {F,LU,TU}
    println(io, "  Restricted Two-body System ", body.name == "" ? "" : string("(",body.name,")"), ":")
    println(io, "")
    println(io, "    μ:  ", body.μ, " ", string(LU^3 / TU^2))
end

"""
Struct for storing `Keplerian` states.
"""
struct KeplerianState{F<:AbstractFloat, LU, TU, AU} <: AbstractUnitfulState{F, LU, TU} where {AU<:Unitful.Units{U, NoDims, nothing} where U}
    e::F
    a::F
    i::F
    Ω::F
    ω::F
    ν::F

    function KeplerianState(e::E, a::A, i::I, Ω::O, ω::W, ν::V; lengthunit = u"km", angularunit = u"rad") where {
        E <: Real, A <: Real, I <: Real, O <: Real, W <: Real, V <: Real
    }
        F = promote_type(E, A, I, O, W, V)
        if !(F <: AbstractFloat)
            @warn "Promoted type ($(string(F)) is not a float: defaulting to Float64."
            F = Float64
        end
        return new{F, lengthunit, u"s", angularunit}(e, a, i, Ω, ω, ν)
    end

    function KeplerianState(e::Real, a::Unitful.Length, i::DimensionlessQuantity, 
                            Ω::DimensionlessQuantity, ω::DimensionlessQuantity, 
                            ν::DimensionlessQuantity)

        F = promote_type(typeof(e), typeof(ustrip(a)), typeof(ustrip(i)), typeof(ustrip(Ω)), typeof(ustrip(ω)), typeof(ustrip(ν)))
        if !(F <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            F = Float64
        end
        lengthunit  = unit(a)
        angularunit = u"rad" ∈ unit.((i, Ω, ω, ν)) ? u"rad" : unit(i)
        return KeplerianState(e, ustrip(lengthunit, a), ustrip(angularunit(i)), 
                              ustrip(angularunit(Ω)), ustrip(angularunit(ω)), ustrip(angularunit(ν));
                              lengthunit = lengthunit, angularunit = angularunit) 
    end

end

"""
Returns the dimmensionless unit associated with the Keplerian state.
"""
angularunit(::K) where K <: KeplerianState = K.parameters[4]

AstrodynamicsCore.coordinateframe(::KeplerianState) = AstrodynamicsCore.Bodycentric

"""
An orbital state within the Restricted Two-body Problem.
"""
struct RestrictedTwoBodyState{C<:AbstractConic, F<:AbstractFloat, T<:Union{CartesianState{F,LU,TU,Bodycentric} where {LU,TU}, KeplerianState{F,LU,TU,AU} where {LU,TU,AU}}} <: AbstractOrbitalState
    state::T
    system::RestrictedTwoBodySystem{F}
    
    function RestrictedTwoBodyState(r, v, body::RestrictedTwoBodySystem{T}, epoch=0) where T <: AbstractFloat
        F = promote_type(eltype(r), eltype(v), T, typeof(epoch))
        if !(F <: AbstractFloat)
            @warn "Promoted type $(string(F)) is not of type float. Defaulting to Float64."
            F = Float64
        end
        state = CartesianState(F.(r), F.(v), epoch, Bodycentric)
        lengthunit = typeof(state).parameters[2]
        timeunit   = typeof(state).parameters[3]
        newbody = RestrictedTwoBodySystem(F(ustrip(lengthunit^3 / timeunit^2, mass_parameter(body))), body.name; lengthunit = lengthunit, timeunit = timeunit)
        return new{conic(eccentricity(state.r, state.v, newbody.μ)), F, typeof(state)}(state, newbody)
    end
    RestrictedTwoBodyState(r, v, μ::Number, epoch=0) = RestrictedTwoBodyState(r, v, RestrictedTwoBodySystem(μ), epoch)

    function RestrictedTwoBodyState(e, a, i, Ω, ω, ν, body::RestrictedTwoBodySystem{T}, epoch=0) where T <: AbstractFloat
        state = KeplerianState(e, a, i, Ω, ω, ν)

        F = typeof(state).parameters[1]
        lengthunit = typeof(state).parameters[2]
        timeunit   = typeof(state).parameters[3]
        newbody = RestrictedTwoBodySystem(F(ustrip(lengthunit^3 / timeunit^2, mass_parameter(body))), body.name; lengthunit = lengthunit, timeunit = timeunit)
        return new{conic(e), F, typeof(state)}(state, newbody)
    end
    RestrictedTwoBodyState(e, a, i, Ω, ω, ν, μ::Number, epoch=0) = RestrictedTwoBodyState(e, a, i, Ω, ω, ν, RestrictedTwoBodySystem(μ), epoch)

end

"""
Alias for `RestrictedTwoBodyState`.
"""
Orbit(r, v, body) = RestrictedTwoBodyState(r, v, body)
Orbit(e, a, i, Ω, ω, ν, body) = RestrictedTwoBodyState(e, a, i, Ω, ω, ν, body)

"""
Custom display for KeplerianState instances.
"""
function Base.show(io::IO, orbit::KeplerianState{F,LU,TU,AU}) where {F,LU,TU,AU}

    println(io, "  Keplerian State:")
    println(io, "")
    println(io, "    e:  ", 
                orbit.e)    
    println(io, "    a:  ", 
                orbit.a, " ", string(LU))
    println(io, "    i:  ", 
                orbit.i, AU == u"rad" ? " " : "", string(AU))
    println(io, "    Ω:  ", 
                orbit.Ω, AU == u"rad" ? " " : "", string(AU))
    println(io, "    ω:  ", 
                orbit.ω, AU == u"rad" ? " " : "", string(AU))
    println(io, "    ν:  ", 
                orbit.ν, AU == u"rad" ? " " : "", string(AU))

end

"""
Custom display for `RestrictedTwoBodyState` instances.
"""
function Base.show(io::IO, orbit::RestrictedTwoBodyState{C, F, T}) where {C, F, T}
    print(io, string(C), " Restricted Two-body Orbit")
    println(" (", string(F), "):")
    println(io, "")
    show(io, orbit.state)
    println(io, "")
    show(io, orbit.system)
end

# Constants

const KeplerianOrbit = RestrictedTwoBodyState{C, F, K} where {C,F,K <: KeplerianState{F, LU, TU, AU} where {LU,TU,AU}}
const CartesianOrbit = RestrictedTwoBodyState{C, F, R} where {C,F,R <: CartesianState{F, LU, TU} where {LU,TU}}

KeplerianOrbit(e, a, i, Ω, ω, ν, system, epoch) = RestrictedTwoBodyState(e, a, i, Ω, ω, ν, system, epoch)
CartesianOrbit(r, v, system, epoch) = RestrictedTwoBodyState(r, v, system, epoch)
KeplerianOrbit(orbit::CartesianOrbit) = KeplerianOrbit(keplerian(position_vector(orbit), velocity_vector(orbit), mass_parameter(orbit.system))..., orbit.system, orbit.state.epoch)
CartesianOrbit(orbit::KeplerianOrbit) = CartesianOrbit(cartesian(eccentricity(orbit), semimajor_axis(orbit), inclination(orbit), RAAN(orbit), argument_of_periapsis(orbit), true_anomoly(orbit), mass_parameter(orbit.system))..., orbit.system, orbit.epoch)