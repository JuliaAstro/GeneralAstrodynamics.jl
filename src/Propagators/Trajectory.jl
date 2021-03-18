# 
# Describes trajectories for all `OrbitalSystem`s.
#

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

Trajectory(traj::Trajectory) = Trajectory(traj.step, traj.t, traj.status)
Base.length(traj::Trajectory) = length(traj.t)
Base.getindex(traj::Trajectory, i) = traj.step[i]
Base.show(io::IO, traj::Trajectory) = println(io, typeof(traj), " with ", length(traj), " steps")
