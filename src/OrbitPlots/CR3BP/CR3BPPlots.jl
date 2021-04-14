#
# Plot Circular Restricted Three-body Problem trajectories
# 

"""
Plot the orbital positions of a CR3BP orbit.
"""
function plotpositions(traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; lengthunit = lengthunit(first(traj).system), kwargs...)

    all_true_or_false(vec) = all(vec) || all(map(v->!v, vec))
    @assert all_true_or_false(map(orbit -> coordinateframe(orbit.state), traj) isa Union{Synodic, Inertial}) "All coordinate frames within a trajectory must be identical for plotting!"
    @assert all_true_or_false(map(orbit -> orbit isa NormalizedCartesianState, traj)) "If one orbit `state` is normalized, all must be!"
    @assert all([isequal(traj[i].system, traj[i+1].system) for i ∈ 1:length(traj)-1]) "All orbits with the trajectory must have the same CR3BP system."

    isnormalized = first(traj).state isa NormalizedCartesianState
    frame        = coordinateframe(first(traj).state)

    pos = transpose(hcat(position_vector.(traj)...))

    if !isnormalized 
        pos = ustrip.(lengthunit, pos)
        LU  = string(lengthunit)
    else
        LU  = string(normalized_length_unit(first(traj).system))
    end

    defaults = (; title     = "$(string(frame)) Positions" * (first(traj).system.name == "" ? "" : " w.r.t. $(first(traj).system.name) Barycenter"),
                  xlabel    = "X ($LU)", 
                  ylabel    = "Y ($LU)", 
                  zlabel    = "Z ($LU)",
                  label     = "Orbit Position",
                  formatter = :plain,
                  grid      = :on,
                  linewidth = 2)

    options = merge(defaults, kwargs)

    plot(pos[:,1], pos[:,2], pos[:,3]; options...)

end

"""
Plot the orbital positions of a CR3BP orbit to the last plot.
"""
function plotpositions!(traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; lengthunit = lengthunit(first(traj).system), kwargs...)

    all_true_or_false(vec) = all(vec) || all(map(v->!v, vec))
    @assert all_true_or_false(map(orbit -> coordinateframe(orbit.state), traj) isa Union{Synodic, Inertial}) "All coordinate frames within a trajectory must be identical for plotting!"
    @assert all_true_or_false(map(orbit -> orbit isa NormalizedCartesianState, traj)) "If one orbit `state` is normalized, all must be!"
    @assert all([isequal(traj[i].system, traj[i+1].system) for i ∈ 1:length(traj)-1]) "All orbits with the trajectory must have the same CR3BP system."

    isnormalized = first(traj).state isa NormalizedCartesianState
    frame        = coordinateframe(first(traj).state)

    pos = transpose(hcat(position_vector.(traj)...))

    if !isnormalized 
        pos = ustrip.(lengthunit, pos)
        LU  = string(lengthunit)
    else
        LU  = string(normalized_length_unit(first(traj).system))
    end

    defaults = (; label = :none)
    options = merge(defaults, kwargs)

    plot!(pos[:,1], pos[:,2], pos[:,3]; options...)

end

"""
Plot the orbital velocities of a CR3BP orbit.
"""
function plotvelocities(traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; velocityunit = velocityunit(first(traj).system), kwargs...)

    all_true_or_false(vec) = all(vec) || all(map(v->!v, vec))
    @assert all_true_or_false(map(orbit -> coordinateframe(orbit.state), traj) isa Union{Synodic, Inertial}) "All coordinate frames within a trajectory must be identical for plotting!"
    @assert all_true_or_false(map(orbit -> orbit isa NormalizedCartesianState, traj)) "If one orbit `state` is normalized, all must be!"
    @assert all([isequal(traj[i].system, traj[i+1].system) for i ∈ 1:length(traj)-1]) "All orbits with the trajectory must have the same CR3BP system."

    isnormalized = first(traj).state isa NormalizedCartesianState
    frame        = coordinateframe(first(traj).state)

    vel = transpose(hcat(velocity_vector.(traj)...))
    
    if !isnormalized 
        vel = ustrip.(velocityunit, vel)
        VU  = string(velocityunit)
    else
        VU  = string(normalized_length_unit(first(traj).system) / normalized_time_unit(first(traj).system))
    end

    defaults = (; title     = "$(string(frame)) Velocities" * (first(traj).system.name == "" ? "" : " w.r.t. $(first(traj).system.name) Barycenter"),
                  xlabel    = "X ($VU)", 
                  ylabel    = "Y ($VU)", 
                  zlabel    = "Z ($VU)",
                  label     = "Orbit Velocity",
                  formatter = :plain,
                  grid      = :on,
                  linewidth = 2)

    options = merge(defaults, kwargs)

    plot(vel[:,1], vel[:,2], vel[:,3]; options...)

end

"""
Plot the orbital velocities of a CR3BP orbit to the last plot.
"""
function plotvelocities!(traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; velocityunit = velocityunit(first(traj).system), kwargs...)

    all_true_or_false(vec) = all(vec) || all(map(v->!v, vec))
    @assert all_true_or_false(map(orbit -> coordinateframe(orbit.state), traj) isa Union{Synodic, Inertial}) "All coordinate frames within a trajectory must be identical for plotting!"
    @assert all_true_or_false(map(orbit -> orbit isa NormalizedCartesianState, traj)) "If one orbit `state` is normalized, all must be!"
    @assert all([isequal(traj[i].system, traj[i+1].system) for i ∈ 1:length(traj)-1]) "All orbits with the trajectory must have the same CR3BP system."

    isnormalized = first(traj).state isa NormalizedCartesianState
    frame        = coordinateframe(first(traj).state)

    vel = transpose(hcat(velocity_vector.(traj)...))
    
    if !isnormalized 
        vel = ustrip.(velocityunit, vel)
        VU  = string(velocityunit)
    else
        VU  = string(normalized_length_unit(first(traj).system) / normalized_time_unit(first(traj).system))
    end

    defaults = (; label = :none)
    options = merge(defaults, kwargs)

    plot(vel[:,1], vel[:,2], vel[:,3]; options...)

end