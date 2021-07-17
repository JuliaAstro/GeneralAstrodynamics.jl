#
# Orbit descriptions.
#

"""
$(TYPEDEF)

A Restricted Two-body Problem system.
"""

"""
$(TYPEDEF)

A supertype for all single-point orbit descriptions.

!!! note
    All `AbstractOrbit` instances __must__ have 
    the following fields: (:epoch, :state, :system).
"""
abstract type AbstractOrbit{F, MU, LU, TU, AU, E, S<:StateVector, P<:ParameterVector} end

"""
$(SIGNATURES)

Returns the epoch (timestamp) for the `AbstractOrbit`.
"""
epoch(orbit::AbstractOrbit) = orbit.epoch

"""
$(SIGNATURES)

Returns the state vector for the `AbstractOrbit`.
"""
state(orbit::AbstractOrbit) = orbit.state

"""
$(SIGNATURES)

Returns the parameter vector for the `AbstractOrbit`.
"""
system(orbit::AbstractOrbit) = orbit.system

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
$(TYPEDEF)

An alias for `Orbit` instances about `R2BP` systems.
"""
const R2BPOrbit = Orbit{F, MU, LU, TU, AU, E, <:Union{CartesianState, KeplerianState}, <:R2BPParameters} where {F, MU, LU, TU, AU, E}

"""
$(TYPEDEF)

An alias for `Orbit` instances about `R2BP` systems.
"""
const CR3BPOrbit = Orbit{F, MU, LU, TU, AU, E, <:CartesianState, <:CR3BPParameters} where {F, MU, LU, TU, AU, E}

"""
$(SIGNATURES)

Outer constructor for `Orbit`s.
"""
function Orbit(state::StateVector, system::ParameterVector, epoch=UTCEpoch(now())) 
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