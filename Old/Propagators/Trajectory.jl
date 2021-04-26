# 
# Describes trajectories for all `AbstractOrbit`s.
#

"""
An alias for a `Vector` of `AbstractOrbits`.
"""
const Trajectory{T} = Vector{T} where T <: AbstractOrbit

Base.show(io::IO, traj::Trajectory{<:CartesianOrbit}) = println(io, "Cartesian Two-body trajectory with ", length(traj), " steps")
Base.show(io::IO, traj::Trajectory{<:KeplerianOrbit}) = println(io, "Keplerian Two-body trajectory with ", length(traj), " steps")
Base.show(io::IO, ::MIME"text/plain", traj::Trajectory) = show(io,traj)
