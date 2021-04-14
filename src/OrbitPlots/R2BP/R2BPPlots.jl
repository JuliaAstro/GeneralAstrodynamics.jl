#
# Plot Restricted Two-body Problem trajectories
# 

"""
Plot the orbital positions of a R2BP orbit.
"""
function plotpositions(traj::Trajectory{<:RestrictedTwoBodyOrbit}; lengthunit = lengthunit(first(traj)), kwargs...)

    pos = transpose(hcat(position_vector.(traj)...))
    pos = ustrip.(lengthunit, pos)
    LU  = string(lengthunit)

    defaults = (; title     = "R2BP Positions" * (first(traj).system.name == "" ? "" : " w.r.t. $(first(traj).system.name)'s Center"),
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
Plot the orbital velocities of a R2BP orbit.
"""
function plotvelocities(traj::Trajectory{<:RestrictedTwoBodyOrbit}; velocityunit = velocityunit(first(traj)), kwargs...)

    vel = transpose(hcat(velocity_vector.(traj)...))
    vel = ustrip.(velocityunit, vel)
    VU  = string(velocityunit)

    defaults = (; title     = "R2BP Velocities" * (first(traj).system.name == "" ? "" : " w.r.t. $(first(traj).system.name)'s Center"),
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