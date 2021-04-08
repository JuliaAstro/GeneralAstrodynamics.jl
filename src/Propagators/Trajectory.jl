# 
# Describes trajectories for all `AbstractOrbitalState`s.
#

"""
A structure for storing trajectories of `TwoBodySystem` orbits,
`RestrictedThreeBodySystem` orbits, and `NBodySystem` orbits.
"""
const Trajectory{T<:AbstractOrbitalState} = Vector{T}

function Trajectory(step::AbstractVector{S}, 
                    system::RestrictedTwoBodySystem,
                    t::AbstractVector{<:Number} = [i for i ∈ 1:length(step)],
                    status::Symbol = :notapplicable) where {S <: Union{CartesianState, KeplerianState}}
    @assert length(step) == length(t) "Time vector and state vectors must have the same length!"
    return RestrictedTwoBodyState.(t, step, (system,))
end

Base.show(io::IO, traj::Trajectory{T}) where {T<:AbstractOrbitalState} = println(io, string(T), " trajectory with ", length(traj), " steps")
