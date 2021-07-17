#
# State vector descriptions.
#

"""
$(TYPEDEF)

A supertype for all state representations in astrodynamics.
"""
abstract type StateVector{F, LU, TU, AU, T} <: ParameterizedLabelledArray{F, 1, T, LArray{F, 1, MVector{6, F}, T}} end

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

Returns the angularunit of the state vector.
"""
velocityunit(::StateVector{F, LU, TU, AU}) where {F, LU, TU, AU} = LU/TU

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
$(TYPEDEF)

A Cartesian state vector with length 6. Internally
uses `MVector` and `LVector` to store data.
Data is accessible via labels, which are
`x, y, z, ẋ, ẏ, ż, r, v`, where `r` 
access `x,y,z` and `v` accesses `ẋ, ẏ, ż`.
"""
mutable struct CartesianState{F, LU, TU, AU} <: StateVector{F, LU, TU, AU, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)} 
    __rawdata::LArray{F, 1, MVector{6, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)}
end

"""
$(SIGNATURES)

Constructs a `CartesianState` from provided position and velocity vectors.
"""
function CartesianState(r::AbstractArray, v::AbstractArray; 
                        lengthunit=(unit(eltype(r)) isa Unitful.Length ? unit(eltype(r)) : u"km"),
                        timeunit=(unit(eltype(v)) isa Unitful.Velocity ? lengthunit / unit(eltype(v)) : u"s"),
                        angularunit=u"°")
    rr = (eltype(r) isa Unitful.Length ? ustrip.(lengthunit, r) : r)
    vv = (eltype(v) isa Unitful.Velocity ? ustrip.(lengthunit / timeunit, v) : v)
    F = promote_type(eltype(rr), eltype(vv))
    F isa AbstractFloat || (F = Float64)
    return CartesianState(vcat(rr,vv); lengthunit = lengthunit, timeunit = timeunit, angularunit = angularunit)
end

"""
$(SIGNATURES)

Constructs a `CartesianState`.
"""
function CartesianState(statevector; lengthunit=u"km", timeunit=u"s", angularunit=u"°")
    F = eltype(statevector)
    F isa AbstractFloat || (F = Float64)

    return CartesianState{F, lengthunit, timeunit, angularunit}(
        LArray{F, 1, MVector{6, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)}(
            statevector
        )
    )
end

"""
$(SIGNATURES)

Returns `x`.
"""
get_x(state::CartesianState) = state[1] * lengthunit(state)

"""
$(SIGNATURES)

Returns `y`.
"""
get_y(state::CartesianState) = state[2] * lengthunit(state)

"""
$(SIGNATURES)

Returns `z`.
"""
get_z(state::CartesianState) = state[3] * lengthunit(state)

"""
$(SIGNATURES)

Returns `ẋ`.
"""
get_ẋ(state::CartesianState) = state[4] * velocityunit(state)

"""
$(SIGNATURES)

Returns `ẏ`.
"""
get_ẏ(state::CartesianState) = state[5] * velocityunit(state)

"""
$(SIGNATURES)

Returns `ż`.
"""
get_ż(state::CartesianState) = state[6] * velocityunit(state)

"""
$(SIGNATURES)

Returns `r`.
"""
get_r(state::CartesianState) = state[1:3] * lengthunit(state)

"""
$(SIGNATURES)

Returns `v`.
"""
get_v(state::CartesianState) = state[4:6] * velocityunit(state)

"""
$(SIGNATURES)

Returns the whole state vector, with units.
"""
statevector(state::CartesianState) = MVector{6}(get_x(state), get_y(state), get_z(state), get_ẋ(state), get_ẏ(state), get_ż(state))

"""
$(SIGNATURES)

Displays a `CartesianState`.
"""
function Base.show(io::IO, state::CartesianState; showfloats=true, space="")
    println(io, space, "Cartesian State:", showfloats ? "($(eltype(state)))" : "", "\n")
    println(io, space, "  Position: [$(state[1]), $(state[2]), $(state[3])] $(lengthunit(state))")
    println(io, space, "  Velocity: [$(state[4]), $(state[5]), $(state[6])] $(lengthunit(state)/timeunit(state))")
end

"""
$(SIGNATURES)

Displays a `CartesianState`.
"""
Base.show(io::IO, ::MIME"text/plain", state::CartesianState; showfloats=true, space="") = show(io, state; showfloats=showfloats, space=space)

"""
$(TYPEDEF)

A Keplerian state vector with length 6. Internally
uses `MVector` and `LVector` to store data.
Data is accessible via labels, which are
`e, a, i, Ω, ω, ν`.
"""
mutable struct KeplerianState{F, LU, TU, AU} <: StateVector{F, LU, TU, AU, (e=1, a=2, i=3, Ω=4, ω=5, ν=6)} 
    __rawdata::LArray{F, 1, MVector{6, F}, (e=1, a=2, i=3, Ω=4, ω=5, ν=6)}

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

Returns `e`.
"""
get_e(state::CartesianState) = state[1]

"""
$(SIGNATURES)

Returns `a`.
"""
get_a(state::CartesianState) = state[2] * lengthunit(state)

"""
$(SIGNATURES)

Returns `i`.
"""
get_i(state::CartesianState) = state[3] * angularunit(state)

"""
$(SIGNATURES)

Returns `Ω`.
"""
get_Ω(state::CartesianState) = state[4] * angularunit(state)

"""
$(SIGNATURES)

Returns `ω.
"""
get_ω(state::CartesianState) = state[5] * angularunit(state)

"""
$(SIGNATURES)

Returns `ν`.
"""
get_ν(state::CartesianState) = state[6] * angularunit(state)

"""
$(SIGNATURES)

Returns the whole state vector, with units.
"""
statevector(state::KeplerianState) = MVector{6}(get_e(state), get_a(state), get_i(state), get_Ω(state), get_ω(state), get_ν(state))

"""
$(SIGNATURES)

Displays a `KeplerianState`.
"""
function Base.show(io::IO, state::KeplerianState; showfloats=true, space="")
    println(io, space, "Keplerian State:", showfloats ? "($(eltype(state)))" : "")
    println(io, space, "  (e=$(state[1]), a=$(state[2]) $(lengthunit(state)), i=$(state[3]) $(angularunit(state)), Ω=$(state[4]) $(angularunit(state)), ω=$(state[5]) $(angularunit(state)), ν=$(state[6]) $(angularunit(state))")
end

"""
$(SIGNATURES)

Displays a `KeplerianState`.
"""
Base.show(io::IO, ::MIME"text/plain", state::KeplerianState; showfloats=true, space="") = show(io, state; showfloats=showfloats, space=space)
