#
# Describes Restricted Two-body orbits through Cartesian coordinates and Orbital Elements.
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
struct RestrictedTwoBodySystem{F, LU, TU} <: AbstractSystem{F, LU, TU}
    μ::F
    name::String

    function RestrictedTwoBodySystem(μ::Real, name=""; lengthunit = u"km", timeunit = u"s")
        T = typeof(μ)
        if !(T <: AbstractFloat)
            @warn "Non-float type $(string(T)) provided. Defaulting to Float64."
            T = Float64
        end
        return new{T, typeof(lengthunit), typeof(timeunit)}(T(μ), name)
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

        lengthunit = (Unitful.FreeUnits{(typeof(M).parameters[1][lengthaxis], ), Unitful.𝐋^3})()^(1//3) 
        timeunit   = (Unitful.FreeUnits{(typeof(M).parameters[1][timeaxis],), Unitful.𝐓^2})()^(-1//2)

        return new{T, typeof(lengthunit), typeof(timeunit)}(T(ustrip(lengthunit^3 / timeunit^2, μ)), name)
    end
end

"""
Conversions between `RestrictedTwoBodySystem` instances.
"""
function Base.convert(::Type{RestrictedTwoBodySystem{F, LU, TU}}, sys::RestrictedTwoBodySystem) where {F, LU, TU}
    return RestrictedTwoBodySystem(ustrip(LU()^3 / TU()^2, mass_parameter(sys)), sys.name; lengthunit = LU(), timeunit = TU())
end

"""
Custom display for `RestrictedTwoBodySystem` instances.
"""
function Base.show(io::IO, body::RestrictedTwoBodySystem{F,LU,TU}) where {F,LU,TU}
    println(io, "  Restricted Two-body System ", body.name == "" ? "" : string("(",body.name,")"), ":")
    println(io, "")
    println(io, "    μ:  ", body.μ, " ", string(LU()^3 / TU()^2))
end

"""
Struct for storing `Keplerian` states.
"""
struct KeplerianState{F, LU, TU, AU} <: AbstractState{F, LU, TU, Inertial} where {AU<:Unitful.Units{U, NoDims, nothing} where U}
    t::F
    e::F
    a::F
    i::F
    Ω::F
    ω::F
    ν::F

    function KeplerianState(e::E, a::A, i::I, Ω::O, ω::W, ν::V, epoch::Real=0; lengthunit = u"km", timeunit = u"s", angularunit = u"rad") where {
        E <: Real, A <: Real, I <: Real, O <: Real, W <: Real, V <: Real
    }
        F = promote_type(E, A, I, O, W, V, typeof(epoch))
        if !(F <: AbstractFloat)
            @warn "Promoted type ($(string(F)) is not a float: defaulting to Float64."
            F = Float64
        end
        return new{F, typeof(lengthunit), typeof(u"s"), typeof(angularunit)}(epoch, e, a, i, Ω, ω, ν)
    end

    function KeplerianState(e::Real, a::Unitful.Length, i::DimensionlessQuantity, 
                            Ω::DimensionlessQuantity, ω::DimensionlessQuantity, 
                            ν::DimensionlessQuantity, epoch::Unitful.Time = 0.0u"s")

        F = promote_type(typeof(e), typeof(ustrip(a)), typeof(ustrip(i)), typeof(ustrip(Ω)), typeof(ustrip(ω)), typeof(ustrip(ν)), typeof(ustrip(epoch)))
        if !(F <: AbstractFloat)
            @warn "Non-float parameters provided. Defaulting to Float64."
            F = Float64
        end
        lengthunit  = unit(a)
        angularunit = u"rad" ∈ unit.((i, Ω, ω, ν)) ? u"rad" : unit(i)
        timeunit = unit(epoch)
        return KeplerianState(e, ustrip(lengthunit, a), ustrip(angularunit, i), 
                              ustrip(angularunit, Ω), ustrip(angularunit, ω), ustrip(angularunit, ν), ustrip(timeunit, epoch);
                              lengthunit = lengthunit, timeunit = timeunit, angularunit = angularunit) 
    end

end

"""
Returns `Unitful` Restricted Two-body orbit eccentricity.
"""
eccentricity(kep::KeplerianState) = kep.e

"""
Returns `Unitful` Restricted Two-body orbit semimajor axis.
"""
semimajor_axis(kep::KeplerianState) = kep.a * lengthunit(kep)

"""
Returns `Unitful` Restricted Two-body orbit inclination.
"""
inclination(kep::KeplerianState) = kep.i * angularunit(kep)

"""
Returns `Unitful` Restricted Two-body orbit right ascension of the ascending node (RAAN).
"""
RAAN(kep::KeplerianState) = kep.Ω * angularunit(kep)

"""
Returns `Unitful` Restricted Two-body orbit argument of periapsis.
"""
argument_of_periapsis(kep::KeplerianState) = kep.ω * angularunit(kep)

"""
Returns `Unitful` Restricted Two-body orbit true anomoly.
"""
true_anomoly(kep::KeplerianState) = kep.ν * angularunit(kep)

"""
Returns the timepoint associated with a `KeplerianState`.
"""
epoch(kep::KeplerianState) = kep.t * timeunit(kep)

"""
Returns the dimmensionless angular unit associated with the Keplerian state.
"""
angularunit(::KeplerianState{F, LU, TU, AU}) where {F, LU, TU, AU} = AU()

"""
An orbital state within the Restricted Two-body Problem.
"""
const R2BPState{F,LU,TU} = Union{CartesianState{F, LU, TU, Inertial}, KeplerianState{F, LU, TU, <:Unitful.DimensionlessUnits}} where {F,LU,TU}

"""
A structure which contains __all__ relevant values for Restricted Two-body Problem orbits.
"""
struct RestrictedTwoBodyOrbit{C<:AbstractConic, F, LU, TU, T<:R2BPState{F,LU,TU}} <: AbstractOrbit{F, LU, TU} 
    state::T
    system::RestrictedTwoBodySystem{F,LU,TU}
    
    function RestrictedTwoBodyOrbit(r, v, body::RestrictedTwoBodySystem{T}, epoch=0) where {T}
        F = promote_type(eltype(ustrip.(r)), eltype(ustrip.(v)), T, typeof(ustrip(epoch)))
        if !(F <: AbstractFloat)
            @warn "Promoted type $(string(F)) is not of type float. Defaulting to Float64."
            F = Float64
        end
        state = CartesianState(r, v, epoch, Inertial)
        LU = lengthunit(state)
        TU   = timeunit(state)
        newbody = RestrictedTwoBodySystem(F(ustrip(LU^3 / TU^2, mass_parameter(body))), body.name; lengthunit = LU, timeunit = TU)
        return new{conic(eccentricity(state.r, state.v, newbody.μ)), F, typeof(LU), typeof(TU), typeof(state)}(state, newbody)
    end
    RestrictedTwoBodyOrbit(r, v, μ::Number, epoch=0) = RestrictedTwoBodyOrbit(r, v, RestrictedTwoBodySystem(μ), epoch)

    function RestrictedTwoBodyOrbit(e, a, i, Ω, ω, ν, body::RestrictedTwoBodySystem{T}, epoch=0) where {T}
        state = KeplerianState(e, a, i, Ω, ω, ν, epoch)
        F = typeof(state).parameters[1]
        LU = lengthunit(state)
        TU   = timeunit(state)
        newbody = RestrictedTwoBodySystem(F(ustrip(LU^3 / TU^2, mass_parameter(body))), body.name; lengthunit = LU, timeunit = TU)
        return new{conic(state.e), F, typeof(LU), typeof(TU), typeof(state)}(state, newbody)
    end
    RestrictedTwoBodyOrbit(e, a, i, Ω, ω, ν, μ::Number, epoch=0) = RestrictedTwoBodyOrbit(e, a, i, Ω, ω, ν, RestrictedTwoBodySystem(μ), epoch)

    function RestrictedTwoBodyOrbit(state::CartesianState, sys::RestrictedTwoBodySystem) 
        F = promote_type(eltype(state), eltype(sys))
        LU = lengthunit(state)
        TU = timeunit(state)
        e = eccentricity(position_vector(state), velocity_vector(state), mass_parameter(sys))
        newstate = convert(CartesianState{F, typeof(LU), typeof(TU)}, state)
        newsys   = convert(RestrictedTwoBodySystem{F, typeof(LU), typeof(TU)}, sys)
        return new{conic(e), F, typeof(LU), typeof(TU), typeof(newstate)}(newstate, newsys)
    end

    function RestrictedTwoBodyOrbit(state::KeplerianState, sys::RestrictedTwoBodySystem) 
        F = promote_type(eltype(state), eltype(sys))
        LU = lengthunit(state)
        TU = timeunit(state)
        AU = angularunit(state)
        e = eccentricity(state)
        newstate = convert(KeplerianState{F, LU, TU, AU}, state)
        newsys   = convert(RestrictedTwoBodySystem{F, LU, TU}, sys)
        return new{conic(e), F, typeof(LU), typeof(TU), typeof(newstate)}(newstate, newsys)
    end

end


"""
Convert between `eltype`, `lengthunit`, and `timeunit` values for `KeplerianStates`.
"""
function Base.convert(::Type{KeplerianState{F,LU,TU,AU}}, kep::KeplerianState) where {F,LU,TU,AU}
    e = F(eccentricity(kep))
    a = F(ustrip(LU(), semimajor_axis(kep)))
    i = F(ustrip(AU(), inclination(kep)))
    Ω = F(ustrip(AU(), RAAN(kep)))
    ω = F(ustrip(AU(), argument_of_periapsis(kep)))
    ν = F(ustrip(AU(), true_anomoly(kep)))
    t = F(ustrip(TU(), epoch(kep)))
    return KeplerianState(e, a, i, Ω, ω, ν, t; lengthunit = LU(), timeunit = TU(), angularunit = AU())
end

# Overrides `Unitful` unit conversions for `AbstractUnitfulStructure` instances.
(u::Unitful.LengthFreeUnits)(state::KeplerianState{F,LU,TU}) where {F,LU,TU} = convert(KeplerianState{F, typeof(u), TU}, state)

# Overrides `Unitful` unit conversions for `AbstractUnitfulStructure` instances.
(u::Unitful.TimeFreeUnits)(state::KeplerianState{F,LU,TU}) where {F,LU,TU} = convert(KeplerianState{F, LU, typeof(u)}, state)


"""
`RestrictedTwoBodyOrbit` instanes are likely the most commonly used feature
of this package. To accomadate this, `Orbit` is an alias for common 
`RestrictedTwoBodyOrbit` constructers.
"""
Orbit(r, v, body, t = 0) = CartesianOrbit(r, v, body, t)

"""
`RestrictedTwoBodyOrbit` instanes are likely the most commonly used feature
of this package. To accomadate this, `Orbit` is an alias for common 
`RestrictedTwoBodyOrbit` constructers.
"""
Orbit(e, a, i, Ω, ω, ν, body, t = 0) = KeplerianOrbit(e, a, i, Ω, ω, ν, body, t)

"""
Print `KeplerianState` instances to `io`.
"""
function Base.show(io::IO, orbit::KeplerianState{F,LU,TU,AU}) where {F,LU,TU,AU}

    println(io, "  Keplerian State:")
    println(io, "")
    println(io, "    e:  ", 
                orbit.e)    
    println(io, "    a:  ", 
                orbit.a, " ", string(LU()))
    println(io, "    i:  ", 
                orbit.i, AU == u"rad" ? " " : "", string(AU()))
    println(io, "    Ω:  ", 
                orbit.Ω, AU == u"rad" ? " " : "", string(AU()))
    println(io, "    ω:  ", 
                orbit.ω, AU == u"rad" ? " " : "", string(AU()))
    println(io, "    ν:  ", 
                orbit.ν, AU == u"rad" ? " " : "", string(AU()))

end

"""
Print `RestrictedTwoBodyOrbit` instances to `io`.
"""
function Base.show(io::IO, orbit::RestrictedTwoBodyOrbit{C, F, T}) where {C, F, T}
    print(io, string(C), " Restricted Two-body Orbit")
    println(" (", string(F), "):")
    println(io, "")
    show(io, orbit.state)
    println(io, "")
    show(io, orbit.system)
end

"""
An alias for `RestrictedTwoBodyOrbit`.
"""
const R2BPOrbit = RestrictedTwoBodyOrbit

"""
An alias for `RestrictedTwoBodySystem`.
"""
const R2BPSystem = RestrictedTwoBodySystem

"""
An alias for `RestrictedTwoBodyOrbit` instances with `KeplerianState` values.
"""
const KeplerianOrbit{C,F,LU,TU} = RestrictedTwoBodyOrbit{C, F, LU, TU, K} where {K <: KeplerianState{F, LU, TU, AU} where AU}

"""
An alias for `RestrictedTwoBodyOrbit` instances with `CartesianState` values.
"""
const CartesianOrbit{C,F,LU,TU} = RestrictedTwoBodyOrbit{C, F, LU, TU, R} where {R <: CartesianState{F, LU, TU, Inertial}}

"""
Alias for a `RestrictedTwoBodyOrbit` constructor with `KeplerianState` values.
"""
KeplerianOrbit(e, a, i, Ω, ω, ν, system, epoch) = RestrictedTwoBodyOrbit(e, a, i, Ω, ω, ν, system, epoch)

"""
Alias for a `RestrictedTwoBodyOrbit` constructor with `CartesianState` values.
"""
CartesianOrbit(r, v, system, epoch) = RestrictedTwoBodyOrbit(r, v, system, epoch)

"""
Alias for a `RestrictedTwoBodyOrbit` constructor with `KeplerianState` values.
"""
KeplerianOrbit(orbit::CartesianOrbit) = KeplerianOrbit(keplerian(position_vector(orbit), velocity_vector(orbit), orbit.system)..., orbit.system, epoch(orbit.state))

"""
Alias for a `RestrictedTwoBodyOrbit` constructor with `CartesianState` values.
"""
CartesianOrbit(orbit::KeplerianOrbit) = CartesianOrbit(cartesian(eccentricity(orbit), semimajor_axis(orbit), inclination(orbit), RAAN(orbit), argument_of_periapsis(orbit), true_anomoly(orbit), mass_parameter(orbit.system))..., orbit.system, epoch(orbit.state))

"""
Alias for a `RestrictedTwoBodyOrbit` constructor with `KeplerianState` values.
"""
KeplerianOrbit(state::KeplerianState, system) = RestrictedTwoBodyOrbit(state, system)

"""
Alias for a `RestrictedTwoBodyOrbit` constructor with `CartesianState` values.
"""
CartesianOrbit(state::CartesianOrbit, system) = RestrictedTwoBodyOrbit(state, system)