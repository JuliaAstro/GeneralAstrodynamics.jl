#
# State vector descriptions.
#

"""
A supertype for all states in astrodynamics.
"""
abstract type AbstractState{F, LU, TU, AU, T} <: ParameterizedLabelledArray{F, 1, T, LArray{F, 1, MVector{6, F}, T}} end

"""
A supertype for all state representations in astrodynamics.
"""
abstract type StateVector{F, LU, TU, AU, T} <: AbstractState{F, LU, TU, AU, T} end

"""
A supertype for all state representations with local linearizations in astrodynamics.
"""
abstract type StateVectorWithSTM{F, LU, TU, AU, T} <: AbstractState{F, LU, TU, AU, T} end

"""
Returns the lengthunit of the state vector.
"""
Base.@pure lengthunit(::AbstractState{F, LU, TU, AU}) where {F, LU, TU, AU} = LU

"""
Returns the timeunit of the state vector.
"""
Base.@pure timeunit(::AbstractState{F, LU, TU, AU}) where {F, LU, TU, AU} = TU

"""
Returns the angularunit of the state vector.
"""
Base.@pure angularunit(::AbstractState{F, LU, TU, AU}) where {F, LU, TU, AU} = AU

"""
Returns the angularunit of the state vector.
"""
velocityunit(::AbstractState{F, LU, TU, AU}) where {F, LU, TU, AU} = LU/TU

"""
The length of any `StateVector` is 6!
"""
Base.@pure Base.length(::StateVector) = 6

"""
The size of any `StateVector` is (6,)!
"""
Base.@pure Base.size(::StateVector) = (6,)

"""
The length of any `StateVector` is 6!
"""
Base.@pure Base.length(::StateVectorWithSTM) = 42

"""
The size of any `StateVector` is (6,)!
"""
Base.@pure Base.size(::StateVectorWithSTM) = (42,)

"""
A Cartesian state vector with length 6. Internally
uses `MVector` and `LVector` to store data.
Data is accessible via labels, which are
`x, y, z, ẋ, ẏ, ż, r, v`, where `r` 
access `x,y,z` and `v` accesses `ẋ, ẏ, ż`.
"""
mutable struct CartesianState{F, LU, TU, AU} <: StateVector{F, LU, TU, AU, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)} 
    __rawdata::LArray{F, 1, MVector{6, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)}


    """
    $(SIGNATURES)
    
    Constructs a `CartesianState` with types specified.
    """
    CartesianState{F, LU, TU, AU}(data) where {F, LU, TU, AU} = new{F,LU,TU,AU}(LArray{F, 1, MVector{6, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6)}(data))
end

"""
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
A Cartesian state vector with local linearization
and length 42. Internally
uses `MVector` and `LVector` to store data.
Data is accessible via labels, which are
`x, y, z, ẋ, ẏ, ż, r, v`, where `r` 
access `x,y,z` and `v` accesses `ẋ, ẏ, ż`.
"""
mutable struct CartesianStateWithSTM{F, LU, TU, AU} <: StateVectorWithSTM{F, LU, TU, AU, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6, Φ=(7:42))} 
    __rawdata::LArray{F, 1, MVector{42, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6, Φ=(7:42))}

    """
    $(SIGNATURES)
    
    Constructs a `CartesianStateWithSTM` with types specified.
    """
    CartesianStateWithSTM{F, LU, TU, AU}(data) where {F, LU, TU, AU} = new{F,LU,TU,AU}(LArray{F, 1, MVector{42, F}, (x = 1, y = 2, z = 3, ẋ = 4, ẏ = 5, ż = 6, r = 1:3, v = 4:6, Φ=(7:42))}(data))
end

"""
Outer constructor for `CartesianStateWithSTM`.
"""
CartesianStateWithSTM(data; lengthunit=u"km", timeunit=u"s", angularunit=u"°") = CartesianStateWithSTM{eltype(data), lengthunit, timeunit, angularunit}(data)

"""
Returns a `CartesianStateWithSTM`, given a `CartesianState`.
"""
CartesianStateWithSTM(cart::CartesianState, stm=Matrix{eltype(cart)}(I(6))) = CartesianStateWithSTM{eltype(cart), lengthunit(cart), timeunit(cart), angularunit(cart)}(cart..., stm...)

"""
All Cartesian state!
"""
const CartesianStateVector = Union{<:CartesianState, <:CartesianStateWithSTM}

"""
Returns `x`.
"""
get_x(state::CartesianStateVector) = state[1]

"""
Returns `y`.
"""
get_y(state::CartesianStateVector) = state[2]

"""
Returns `z`.
"""
get_z(state::CartesianStateVector) = state[3]

"""
Returns `ẋ`.
"""
get_ẋ(state::CartesianStateVector) = state[4]

"""
Returns `ẏ`.
"""
get_ẏ(state::CartesianStateVector) = state[5]

"""
Returns `ż`.
"""
get_ż(state::CartesianStateVector) = state[6]

"""
Returns `r`.
"""
get_r(state::CartesianStateVector) = state[1:3]

"""
Returns `v`.
"""
get_v(state::CartesianStateVector) = state[4:6]

"""
Returns `ϕ`.
"""
get_ϕ(state::CartesianStateWithSTM) = state[7:42]

"""
Returns the state transition matrix.
"""
get_stm(state::CartesianStateWithSTM) = MMatrix{6,6}(state[7:42])

"""
Returns the whole state vector, without units.
"""
statevector(state::CartesianState) = MVector{6}(get_x(state), get_y(state), get_z(state), get_ẋ(state), get_ẏ(state), get_ż(state))

"""
Returns the whole state vector, without units.
"""
statevector(state::CartesianStateWithSTM) = MVector{42}(get_x(state), get_y(state), get_z(state), get_ẋ(state), get_ẏ(state), get_ż(state), get_ϕ(state)...)

"""
Displays a `CartesianState`.
"""
function Base.show(io::IO, state::CartesianState; showfloats=true, space="")
    LU = isnormalized(state) ? one(eltype(state)) : lengthunit(state)
    VU = isnormalized(state) ? one(eltype(state)) : velocityunit(state)
    println(io, space, "Cartesian State", showfloats ? " with eltype $(eltype(state))" : "", "\n")
    println(io, space, "  ", "x = $(get_x(state) * LU)")
    println(io, space, "  ", "y = $(get_y(state) * LU)")
    println(io, space, "  ", "z = $(get_z(state) * LU)")
    println(io, space, "  ", "ẋ = $(get_ẋ(state) * VU)")
    println(io, space, "  ", "ẏ = $(get_ẏ(state) * VU)")
    println(io, space, "  ", "ż = $(get_ż(state) * VU)")
end
Base.show(io::IO, ::MIME"text/plain", state::CartesianState; showfloats=true, space="") = show(io, state; showfloats=showfloats, space=space)

"""
Displays a `CartesianStateWithSTM`.
"""
function Base.show(io::IO, state::CartesianStateWithSTM; showfloats=true, space="")
    LU = isnormalized(state) ? one(eltype(state)) : lengthunit(state)
    VU = isnormalized(state) ? one(eltype(state)) : velocityunit(state)
    println(io, space, "Cartesian State with local linearization", showfloats ? " and eltype $(eltype(state))" : "", "\n")
    println(io, space, "  ", "x = $(get_x(state) * LU)")
    println(io, space, "  ", "y = $(get_y(state) * LU)")
    println(io, space, "  ", "z = $(get_z(state) * LU)")
    println(io, space, "  ", "ẋ = $(get_ẋ(state) * VU)")
    println(io, space, "  ", "ẏ = $(get_ẏ(state) * VU)")
    println(io, space, "  ", "ż = $(get_ż(state) * VU)")
end
Base.show(io::IO, ::MIME"text/plain", state::CartesianStateWithSTM; showfloats=true, space="") = show(io, state; showfloats=showfloats, space=space)


"""
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
Constructs a `KeplerianState` from any `AbstractVector`.
"""
KeplerianState{F, LU, TU, AU}(statevector::AbstractVector{<:Real}) where {F, LU, TU, AU} = KeplerianState(
    convert(LArray{F, 1, MVector{6, F}, (e=1, a=2, i=3, Ω=4, ω=5, ν=6)}, statevector)
)

"""
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
Returns `e`.
"""
get_e(state::KeplerianState) = state[1]

"""
Returns `a`.
"""
get_a(state::KeplerianState) = state[2]

"""
Returns `i`.
"""
get_i(state::KeplerianState) = state[3]

"""
Returns `Ω`.
"""
get_Ω(state::KeplerianState) = state[4]

"""
Returns `ω.
"""
get_ω(state::KeplerianState) = state[5]

"""
Returns `ν`.
"""
get_ν(state::KeplerianState) = state[6]

"""
Returns the whole state vector, without units.
"""
statevector(state::KeplerianState) = MVector{6}(get_e(state), get_a(state), get_i(state), get_Ω(state), get_ω(state), get_ν(state))

"""
Displays a `KeplerianState`.
"""
function Base.show(io::IO, state::KeplerianState; showfloats=true, space="")
    println(io, space, "Keplerian State", showfloats ? " with eltype $(eltype(state))" : "", "\n")
    println(io, space, "  ", "e = $(get_e(state))")
    println(io, space, "  ", "a = $(get_a(state) * lengthunit(state))")
    println(io, space, "  ", "Ω = $(get_Ω(state) * angularunit(state))")
    println(io, space, "  ", "ω = $(get_ω(state) * angularunit(state))")
    println(io, space, "  ", "ν = $(get_ν(state) * angularunit(state))")
end
Base.show(io::IO, ::MIME"text/plain", state::KeplerianState; showfloats=true, space="") = show(io, state; showfloats=showfloats, space=space)

"""
Returns `true` if the `StateVector` has `lengthunit` and `timeunit`
parameters of some type `T <: AbstractQuantity`. This is
intended to be used for normalizing `CartesianState` vectors.
"""
isnormalized(state::StateVector) = lengthunit(state) isa Unitful.AbstractQuantity && timeunit(state) isa Unitful.AbstractQuantity

"""
Converts types and units for a `CartesianState`.
"""
function Base.convert(::Type{CartesianState{F, LU, TU, AU}}, state::CartesianState) where {F, LU, TU, AU}
    r = States.get_r(state) .* lengthunit(state) ./ LU .|> upreferred .|> F
    v = States.get_v(state) .* velocityunit(state) ./ (LU/TU) .|> upreferred .|> F
    return CartesianState(vcat(r,v); lengthunit=LU, timeunit=TU, angularunit=AU)
end

"""
Converts types and units for a `CartesianState`.
"""
function Base.convert(::Type{KeplerianState{F, LU, TU, AU}}, state::KeplerianState) where {F, LU, TU, AU}
    e = get_e(state) |> F
    a = ustrip(LU, get_a(state) * lengthunit(state))  |> F
    i = ustrip(AU, get_i(state) * angularunit(state)) |> F
    Ω = ustrip(AU, get_Ω(state) * angularunit(state)) |> F
    ω = ustrip(AU, get_ω(state) * angularunit(state)) |> F
    ν = ustrip(AU, get_ν(state) * angularunit(state)) |> F
    return KeplerianState(MVector{6}(e, a, i, Ω, ω, ν); lengthunit=LU, timeunit=TU, angularunit=AU)
end