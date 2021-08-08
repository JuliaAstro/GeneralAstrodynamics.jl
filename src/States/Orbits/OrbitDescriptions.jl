#
# Orbit descriptions.
#

"""
A supertype for all single-point orbit descriptions.
Parameterized by coordinate frame, floating point type,
mass unit, lenth unit, time unit, epoch type, state type,
and parameter type (in order).
"""
abstract type AbstractOrbit{FR, F, MU, LU, TU, AU, E, S<:AbstractState{F, LU, TU, AU}, P<:ParameterVector{F, MU, LU, TU, AU}} end

"""
An orbit, described by a `StateVector`, with parameters described by a `ParameterVector`.
"""
struct Orbit{FR, F, MU, LU, TU, AU, E, S, P} <: AbstractOrbit{FR, F, MU, LU, TU, AU, E, S, P}
    epoch::E
    state::S
    system::P
end

"""
Returns a default frame.
"""
function defaultframe(system::ParameterVector)
    if system isa R2BPParameters
        return CoordinateFrames.BodycentricInertial
    elseif system isa CR3BPParameters
        return CoordinateFrames.BarycentricRotating
    else
        return Inertial
    end
end

"""
Outer constructor for `Orbit`s.
"""
function Orbit(state::AbstractState, system::ParameterVector, epoch=TAIEpoch(now()); frame = defaultframe(system)) 

    S = eval(Base.typename(typeof(state)).name)
    P = eval(Base.typename(typeof(system)).name)

    if state isa KeplerianState && !(system isa R2BPParameters)
        throw(ArgumentError("Orbits with $P must have CartesianState state vectors."))
    end

    FR = frame
    F  = promote_type(eltype(state), eltype(system))
    MU = massunit(system)

    if system isa CR3BPParameters
        LU = lengthunit(system)
        TU = timeunit(system)
    else
        LU = lengthunit(state)
        TU = timeunit(state)
    end

    AU = angularunit(state)
    E  = typeof(epoch)

    return Orbit{FR, F, MU, LU, TU, AU, E, S{F, LU, TU, AU}, P{F, MU, LU, TU, AU}}(
        epoch,
        convert(S{F, LU, TU, AU}, state),
        convert(P{F, MU, LU, TU, AU}, system)
    )
end

"""
Shows all `Orbit` instances.
"""
function Base.show(io::IO, orbit::Orbit)
    if system(orbit) isa R2BPParameters
        println(io, conic(orbit), " Restricted Two-body Orbit with eltype $(typeof(orbit).parameters[2])")
    elseif system(orbit) isa CR3BPParameters
        println(io, "Circular Restricted Three-body Orbit with eltype $(typeof(orbit).parameters[2])")
    end

    println(io, "")
    Base.show(io, state(orbit); showfloats=false, space="  ")
    println(io, "")
    Base.show(io, system(orbit); showfloats=false, space="  ")
    
end
Base.show(io::IO, ::MIME"text/plain", orbit::Orbit) = show(io, orbit)

"""
Returns the epoch (timestamp) for the `Orbit`.
"""
epoch(orbit::Orbit) = orbit.epoch

"""
Returns the state vector for the `Orbit`.
"""
state(orbit::Orbit) = orbit.state

"""
Returns the parameter vector for the `Orbit`.
"""
system(orbit::Orbit) = orbit.system

"""
Returns the `OrbitalFrame` for the `Orbit`.
"""
Base.@pure frame(::Orbit{FR}) where FR = FR

"""
An alias for `Orbit` instances about `R2BP` systems.
"""
const R2BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:Union{CartesianStateVector, KeplerianState}, <:R2BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about `R2BP` systems with `KeplerianState` descriptions.
"""
const KeplerianR2BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:KeplerianState, <:R2BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about `R2BP` systems with `CartesianState` descriptions.
"""
const CartesianR2BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:CartesianStateVector, <:R2BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about `CR3BP` systems.
"""
const CR3BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:CartesianStateVector, <:CR3BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about any systems with `CartesianState` descriptions.
"""
const CartesianOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:CartesianStateVector, <:ParameterVector} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about any systems with `CartesianStateWithSTM` descriptions.
"""
const CartesianOrbitWithSTM = Orbit{FR, F, MU, LU, TU, AU, E, <:CartesianStateWithSTM, <:ParameterVector} where {FR, F, MU, LU, TU, AU, E}