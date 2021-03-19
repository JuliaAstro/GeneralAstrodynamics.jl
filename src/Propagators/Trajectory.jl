# 
# Describes trajectories for all `OrbitalSystem`s.
#

"""
A structure for storing trajectories of `TwoBodySystem` orbits,
`RestrictedThreeBodySystem` orbits, and `NBodySystem` orbits.
"""
struct Trajectory{T<:OrbitalSystem} <: AbstractTrajectory
    t::Vector{<:Number}
    step::Vector{T}
    status::Symbol

    function Trajectory(step::AbstractVector{T}, 
                        t::AbstractVector{<:Number} = [i for i ∈ 1:length(step)],
                        status::Symbol = :notapplicable) where T <: OrbitalSystem
        @assert length(step) == length(t) "Time vector and state vectors must have the same length!"
        return new{T}(Vector(t), Vector(step), status)
    end
end

"""
Copy constructor for `Trajectory` instances.
"""
Trajectory(traj::Trajectory) = Trajectory(traj.step, traj.t, traj.status)

"""
The `length` of a trajectory is the number of steps in the trajectory.
"""
Base.length(traj::Trajectory) = length(traj.t)

"""
The _n-th_ `index` of a trajectory is the _n-th_ step of the trajectory.
"""
Base.getindex(traj::Trajectory, i) = traj.step[i]
Base.show(io::IO, traj::Trajectory) = println(io, typeof(traj), " with ", length(traj), " steps")
