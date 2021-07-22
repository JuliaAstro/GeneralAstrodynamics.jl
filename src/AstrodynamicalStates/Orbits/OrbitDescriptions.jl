#
# Orbit descriptions.
#

"""
A supertype for all single-point orbit descriptions.
Parameterized by coordinate frame, floating point type,
mass unit, lenth unit, time unit, epoch type, state type,
and parameter type (in order).
"""
abstract type AbstractOrbit{FR, F, MU, LU, TU, AU, E, S<:StateVector, P<:ParameterVector} end

"""
An orbit, described by a `StateVector`, with parameters described by a `ParameterVector`.
"""
struct Orbit{FR, F, MU, LU, TU, AU, E, S, P} <: AbstractOrbit{FR, F, MU, LU, TU, AU, E, S, P}
    epoch::E
    state::S
    system::P
end

"""
Outer constructor for `Orbit`s.
"""
function Orbit(state::StateVector, system::ParameterVector, epoch=UTCEpoch(now()); frame=Inertial) 
    FR = frame
    F  = promote_type(eltype(state), eltype(system))
    MU = massunit(system)
    LU = lengthunit(state)
    TU = timeunit(state)
    AU = angularunit(state)
    E  = typeof(epoch)
    S  = typeof(state)
    P  = typeof(system)

    return Orbit{FR, F,MU,LU,TU,AU,E,S,P}(
        epoch,
        state,
        system
    )
end

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
An alias for `Orbit` instances about `R2BP` systems.
"""
const R2BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:Union{CartesianState, KeplerianState}, <:R2BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about `R2BP` systems with `KeplerianState` descriptions.
"""
const KeplerianR2BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:KeplerianState, <:R2BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about `R2BP` systems with `CartesianState` descriptions.
"""
const CartesianR2BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:CartesianState, <:R2BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about `CR3BP` systems.
"""
const CR3BPOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:CartesianState, <:CR3BPParameters} where {FR, F, MU, LU, TU, AU, E}

"""
An alias for `Orbit` instances about any systems with `CartesianState` descriptions.
"""
const CartesianOrbit = Orbit{FR, F, MU, LU, TU, AU, E, <:CartesianState, <:ParameterVector} where {FR, F, MU, LU, TU, AU, E}
