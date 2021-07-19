#
# Orbit descriptions.
#

"""
$(TYPEDEF)

A supertype for all single-point orbit descriptions.
"""
abstract type AbstractOrbit{F, MU, LU, TU, AU, E, S<:StateVector, P<:ParameterVector} end

"""
$(TYPEDEF)

An orbit, described by a `StateVector`, with parameters described by a `ParameterVector`.
"""
struct Orbit{F, MU, LU, TU, AU, E, S, P} <: AbstractOrbit{F, MU, LU, TU, AU, E, S, P}
    epoch::E
    state::S
    system::P
end

"""
$(SIGNATURES)

Outer constructor for `Orbit`s.
"""
function Orbit(state::StateVector, system::ParameterVector, epoch=AstroTime.UTCEpoch(now())) 
    F  = promote_type(eltype(state), eltype(system))
    MU = massunit(system)
    LU = lengthunit(state)
    TU = timeunit(state)
    AU = angularunit(state)
    E  = typeof(epoch)
    S  = typeof(state)
    P  = typeof(system)

    return Orbit{F,MU,LU,TU,AU,E,S,P}(
        epoch,
        state,
        system
    )
end

"""
$(SIGNATURES)

Returns the epoch (timestamp) for the `Orbit`.
"""
epoch(orbit::Orbit) = orbit.epoch

"""
$(SIGNATURES)

Returns the state vector for the `Orbit`.
"""
state(orbit::Orbit) = orbit.state

"""
$(SIGNATURES)

Returns the parameter vector for the `Orbit`.
"""
system(orbit::Orbit) = orbit.system

"""
$(TYPEDEF)

An alias for `Orbit` instances about `R2BP` systems.
"""
const R2BPOrbit = Orbit{F, MU, LU, TU, AU, E, <:Union{CartesianState, KeplerianState}, <:R2BPParameters} where {F, MU, LU, TU, AU, E}

"""
$(TYPEDEF)

An alias for `Orbit` instances about `CR3BP` systems.
"""
const CR3BPOrbit = Orbit{F, MU, LU, TU, AU, E, <:CartesianState, <:CR3BPParameters} where {F, MU, LU, TU, AU, E}