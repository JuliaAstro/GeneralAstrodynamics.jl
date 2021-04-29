#
# Plot Circular Restricted Three-body Problem trajectories
# 

"""
Plot the orbital positions of a CR3BP orbit.
"""
function plotpositions(traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; lengthunit = lengthunit(first(traj).system), exclude_z = false, kwargs...)

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
                  linewidth = 2,
                  dpi       = 150)

    options = merge(defaults, kwargs)

    if !exclude_z
        return plot(pos[:,1], pos[:,2], pos[:,3]; options...)
    else 
        return plot(pos[:,1], pos[:,2]; options...)
    end
end

"""
Plot the orbital positions of a CR3BP orbit to the last plot.
"""
function plotpositions!(fig, traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; lengthunit = lengthunit(first(traj).system), exclude_z = false, kwargs...)

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

    if !exclude_z
        return plot!(fig, pos[:,1], pos[:,2], pos[:,3]; options...)
    else 
        return plot!(fig, pos[:,1], pos[:,2]; options...)
    end
end

"""
Plot the orbital velocities of a CR3BP orbit.
"""
function plotvelocities(traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; velocityunit = velocityunit(first(traj).system), exclude_z = false, kwargs...)

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
                  linewidth = 2,
                  dpi       = 150)
                                    
    options = merge(defaults, kwargs)

    if !exclude_z
        return plot(vel[:,1], vel[:,2], vel[:,3]; options...)
    else 
        return plot(vel[:,1], vel[:,2]; vel...)
    end
end

"""
Plot the orbital velocities of a CR3BP orbit to the last plot.
"""
function plotvelocities!(fig, traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}; velocityunit = velocityunit(first(traj).system), exclude_z = false, kwargs...)

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

    if !exclude_z
        return plot!(fig, vel[:,1], vel[:,2], vel[:,3]; options...)
    else 
        return plot!(fig, vel[:,1], vel[:,2]; vel...)
    end
end

"""
Plot the nondimensional positions of 
CR3BP celestial bodies.
"""
function plotbodies(system::CircularRestrictedThreeBodySystem; normalize = true, exclude_z = true, kwargs...)

    μ = normalized_mass_parameter(system)

    if normalize
        r₁ = [[-μ],  [0], [0]]
        r₂ = [[1-μ], [0], [0]]
        LU = string(normalized_length_unit(system))
    else
        r₁ = [[ustrip(length_unit(system), -μ * normalized_length_unit(system))],  [0], [0]]
        r₂ = [[ustrip(length_unit(system), 1-μ * normalized_length_unit(system))], [0], [0]]
        LU = string(lengthunit(system))
    end

    defaults = (; title = "CR3BP Bodies", xlabel = "X ($LU)", ylabel = "Y ($LU)", zlabel = "Z ($LU)")
    options  = merge(defaults, kwargs)

    if exclude_z
        fig = scatter(r₁[1:2]...; markersize=10, markercolor=:lightblue, label="Body 1")
        scatter!(fig, r₂[1:2]...; markersize=6, markercolor=:gray, label="Body 2")
        plot!(fig; options...)
        return fig
    else
        fig = scatter(r₁...; markersize=10, markercolor=:lightblue, label="Body 1")
        scatter!(fig, r₂...; markersize=6, markercolor=:gray, label="Body 2")
        plot!(fig; options...)
        return fig
    end

end


"""
Plot the nondimensional positions of 
CR3BP celestial bodies to the last figure.
"""
function plotbodies!(fig, system::CircularRestrictedThreeBodySystem; normalize = true, exclude_z = true, kwargs...)

    μ = normalized_mass_parameter(system)

    if normalize
        r₁ = [[-μ],  [0], [0]]
        r₂ = [[1-μ], [0], [0]]
        LU = string(normalized_length_unit(system))
    else
        r₁ = [[ustrip(length_unit(system), -μ * normalized_length_unit(system))],  [0], [0]]
        r₂ = [[ustrip(length_unit(system), 1-μ * normalized_length_unit(system))], [0], [0]]
        LU = string(lengthunit(system))
    end

    defaults = (; title = "CR3BP Bodies", xlabel = "X ($LU)", ylabel = "Y ($LU)", zlabel = "Z ($LU)")
    options  = merge(defaults, kwargs)

    if exclude_z
        scatter!(fig, r₁[1:2]...; markersize=10, markercolor=:lightblue, label="Body 1")
        scatter!(fig, r₂[1:2]...; markersize=6, markercolor=:gray, label="Body 2")
        plot!(fig; options...)
        return fig
    else
        scatter!(fig, r₁...; markersize=10, markercolor=:lightblue, label="Body 1")
        scatter!(fig, r₂...; markersize=6, markercolor=:gray, label="Body 2")
        plot!(fig; options...)
        return fig
    end

end