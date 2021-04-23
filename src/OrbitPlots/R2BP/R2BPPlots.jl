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

"""
Plot the positions of an orbit.
"""
function plotpositions(pos::AbstractVector{V}; lengthunit = u"km", exclude_z = false, kwargs...) where V <: AbstractVector{<:Real}
    
    pos = transpose(hcat(pos...))
    LU  = string(lengthunit)

    defaults = (; title     = "Orbit Positions",
                  xlabel    = "X ($LU)", 
                  ylabel    = "Y ($LU)", 
                  zlabel    = "Z ($LU)",
                  label     = "Orbit Position",
                  formatter = :plain,
                  grid      = :on,
                  linewidth = 2)

    options = merge(defaults, kwargs)

    if !exclude_z
        plot(pos[:,1], pos[:,2], pos[:,3]; options...)
    else 
        plot(pos[:,1], pos[:,2]; options...)
    end
end

"""
Plot the positions of an orbit.
"""
function plotpositions(pos::AbstractVector{V}; lengthunit = unit(pos[1][1]), exclude_z = false, kwargs...) where V <: AbstractVector{<:Unitful.Length}
    return plotpositions(ustrip.(lengthunit, pos); lengthunit = lengthunit, exclude_z = exclude_z, kwargs...)
end

"""
Plot the positions of an orbit.
"""
function plotpositions!(fig, pos::AbstractVector{V}; lengthunit = u"km", exclude_z = false, kwargs...) where V <: AbstractVector{<:Real}
    
    pos = transpose(hcat(pos...))
    LU  = string(lengthunit)

    defaults = (; title     = "Orbit Positions",
                  xlabel    = "X ($LU)", 
                  ylabel    = "Y ($LU)", 
                  zlabel    = "Z ($LU)",
                  label     = "Orbit Position",
                  formatter = :plain,
                  grid      = :on,
                  linewidth = 2)

    options = merge(defaults, kwargs)

    if !exclude_z
        plot!(fig, pos[:,1], pos[:,2], pos[:,3]; options...)
    else 
        plot!(fig, pos[:,1], pos[:,2]; options...)
    end
end

"""
Plot the positions of an orbit.
"""
function plotpositions!(fig, pos::AbstractVector{V}; lengthunit = unit(pos[1][1]), exclude_z = false, kwargs...) where V <: AbstractVector{<:Unitful.Length}
    return plotpositions!(fig, ustrip.(lengthunit(pos[1][1]), pos); lengthunit = lengthunit, exclude_z = exclude_z, kwargs...)
end

"""
Plot the positions of an orbit.
"""
function plotpositions(pos::AbstractMatrix{<:Real}; lengthunit = u"km", exclude_z = false, kwargs...) 
    
    size(pos,2) == 3 || (pos = transpose(pos))
    LU  = string(lengthunit)

    defaults = (; title     = "Orbit Positions",
                  xlabel    = "X ($LU)", 
                  ylabel    = "Y ($LU)", 
                  zlabel    = "Z ($LU)",
                  label     = "Orbit Position",
                  formatter = :plain,
                  grid      = :on,
                  linewidth = 2)

    options = merge(defaults, kwargs)

    if !exclude_z
        plot(pos[:,1], pos[:,2], pos[:,3]; options...)
    else 
        plot(pos[:,1], pos[:,2]; options...)
    end
end

"""
Plot the positions of an orbit.
"""
function plotpositions(pos::AbstractMatrix{<:Unitful.Length}; lengthunit = unit(pos[1]), exclude_z = false, kwargs...) 
    return plotpositions(ustrip.(lengthunit(pos[1]), pos); lengthunit = lengthunit, exclude_z = exclude_z, kwargs...)
end

"""
Plot the positions of an orbit.
"""
function plotpositions!(fig, pos::AbstractMatrix{<:Real}; lengthunit = u"km", exclude_z = false, kwargs...) 
    
    size(pos,2) == 3 || (pos = transpose(pos))
    LU  = string(lengthunit)

    defaults = (; title     = "Orbit Positions",
                  xlabel    = "X ($LU)", 
                  ylabel    = "Y ($LU)", 
                  zlabel    = "Z ($LU)",
                  label     = "Orbit Position",
                  formatter = :plain,
                  grid      = :on,
                  linewidth = 2)

    options = merge(defaults, kwargs)

    if !exclude_z
        plot!(fig, pos[:,1], pos[:,2], pos[:,3]; options...)
    else 
        plot!(fig, pos[:,1], pos[:,2]; options...)
    end
end

"""
Plot the positions of an orbit.
"""
function plotpositions!(fig, pos::AbstractMatrix{<:Unitful.Length}; lengthunit = unit(pos[1]), exclude_z = false, kwargs...) 
    return plotpositions!(fig, ustrip.(lengthunit(pos[1]), pos); lengthunit = lengthunit, exclude_z = exclude_z, kwargs...)
end