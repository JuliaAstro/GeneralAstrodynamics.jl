#
# Types for orbital trajectories
#

"""
An alias for some abstract `ODESolution` 
with `States` state vector and parameter 
vector types.
"""
const AbstractOrbitalODESolution = SciMLBase.AbstractODESolution{T,N,<:AbstractVector{U} where U<:States.AbstractState} where {T,N}

"""
An wrapper for a `SciMLBase.ODESolution` with a `GeneralAstrodynamics.States.AbstractState` 
state vector type. This represents an object's `Trajectory` in space!
"""
struct Trajectory{FR, E, S<:AbstractOrbitalODESolution}
    epoch::E
    solution::S
end

"""
Returns the start `epoch`. Typically,
This is a type defined in `AstroTime.Epochs`.
"""
initialepoch(traj::Trajectory) = traj.epoch

"""
Returns the `solution` for the `Trajectory`. Typically, 
this is a `DifferentialEquations.ODESolution`.
"""
solution(traj::Trajectory) = traj.solution

"""
A wrapper for (::ODESolution-like)(args...). Returns a state 
of type `initialstate)` at time `t` past `epoch`.
"""
(traj::Trajectory)(t; continuity=:left) = typeof(initialstate(traj))(traj.solution(t, Val{0}, nothing, continuity))

"""
Returns the initial condition associated with the `Trajectory`.
"""
initialstate(traj::Trajectory) = solution(traj).prob.u0

"""
Returns the system associated with the `Trajectory`.
"""
States.system(traj::Trajectory) = solution(traj).prob.p

"""
Returns the `eltype` of the `Trajectory`.
"""
Base.eltype(traj::Trajectory) = eltype(solution(traj))

"""
The `length` of a `Trajectory`.
"""
Base.length(traj::Trajectory) = length(solution(traj))

"""
Calls the underlying `solution`'s `getindex` function
to return the `CartesianState` of the `Trajectory`
at time `t` past the `initialepoch`.
"""
Base.getindex(traj::Trajectory, args...) = getindex(solution(traj), args...)

"""
The `size` of a `Trajectory`.
"""
Base.size(traj::Trajectory) = size(solution(traj))

"""
Returns the `OrbitalFrame` of the `Trajectory`.
"""
Base.@pure States.frame(::Trajectory{FR}) where FR = FR

"""
Returns the length unit for the `Trajectory`.
"""
States.lengthunit(traj::Trajectory) = lengthunit(initialstate(traj))

"""
Returns the time unit for the `Trajectory`.
"""
States.timeunit(traj::Trajectory) = timeunit(initialstate(traj))

"""
Returns the angular unit for the `Trajectory`.
"""
States.angularunit(traj::Trajectory) = angularunit(initialstate(traj))

"""
Returns the velocity unit for the `Trajectory`.
"""
States.velocityunit(traj::Trajectory) = velocityunit(initialstate(traj))

"""
Returns the mass unit for the `Trajectory`.
"""
States.massunit(traj::Trajectory) = massunit(system(traj))

"""
Returns the mass parameter unit for the `Trajectory`.
"""
States.massparamunit(traj::Trajectory) = massparamunit(system(traj))

"""
Returns the position of the `Trajectory` at `t` `timeunit`'s from
the `initialstate`'s `epoch`.
"""
Base.position(traj::Trajectory, t) = t isa Unitful.Time ? position(traj, t / timeunit(traj)) : traj(t, Val{0}; idxs=1:3, continuity=:left) * lengthunit(traj)

"""
Returns the scalar distance of the `Trajectory` at `t` `timeunit`'s from
the `initialstate`'s `epoch`.
"""
States.distance(traj::Trajectory, t) = t isa Unitful.Time ? distance(traj, t / timeunit(traj)) : norm(position(traj, t))

"""
Returns the position of the `Trajectory` at `t` `timeunit`'s from
the `initialstate`'s `epoch`.
"""
States.velocity(traj::Trajectory, t) = t isa Unitful.Time ? velocity(traj, t / timeunit(traj)) : traj(t, Val{0}; idxs=4:6, continuity=:left) * velocityunit(traj)

"""
Returns the scalar distance of the `Trajectory` at `t` `timeunit`'s from
the `initialstate`'s `epoch`.
"""
States.speed(traj::Trajectory, t) = t isa Unitful.Time ? speed(traj, t / timeunit(traj)) : norm(velocity(traj, t))

"""
Returns a `state` of type `typeof(initialstate)`
at time `t`.
"""
States.state(traj::Trajectory, t) = t isa Unitful.Time ? state(traj, t / timeunit(traj)) : typeof(initialstate(traj))(solution(traj)(t))

"""
Returns an `epoch` at time `t` past `initialepoch`.
"""
States.epoch(traj::Trajectory, t) = t isa Unitful.Time ? initialepoch(traj) + UnitfulToAstroTime(t) : epoch(traj, t * timeunit(traj))

"""
Converts `Unitful` types to `AstroTime` types.
Throws an `ArgumentError` if the unit of quantity 
`t` is not a second, minute, hour, day, or year.
"""
function UnitfulToAstroTime(t::Unitful.Time)
    un  = Unitful.unit(t)
    val = ustrip(un, t)

    if un == u"s"
        return val * AstroTime.seconds
    elseif un == u"minute"
        return val * Astrotime.minutes
    elseif un == u"hr"
        return val * AstroTime.hours
    elseif un == u"d"
        return val * AstroTime.days
    elseif un == u"yr"
        return val * AstroTime.years
    else
        throw(ArgumentError("Invalid time unit!"))
    end
end

"""
Returns an `Orbit` at time `t`.
"""
States.Orbit(traj::Trajectory, t) = Orbit(state(traj, t), system(traj), epoch(traj, t))

"""
Show a `Trajectory`.
"""
function Base.show(io::IO, traj::Trajectory)
    println(io, "Trajectory with $(length(traj)) timesteps and eltype $(eltype(solution(traj)))\n")
end