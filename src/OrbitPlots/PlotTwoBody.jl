#
# Plot orbits
#
# Referencing:
# [1] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

"""
Plots every timestep in `sols` in `3D` space. All keyward 
arguments are passed directly to `Plots.jl`.
"""
function orbitplot(sols::Trajectory{<:RestrictedTwoBodyOrbit}, frame=:Cartesian; kwargs...)
   
    # Provided frame can be :Cartesian, or :Perifocal
    if frame == :Perifocal
        return plot2d(sols; kwargs...)
    elseif frame == :Cartesian
        return plot3d(sols; kwargs...)
    else 
        throw(ArgumentError("`frame` must be set to `:Cartesian` or `:Perifocal`"))
    end

end

function plot2d(sols::Trajectory{<:RestrictedTwoBodyOrbit}; kwargs...)

    # Set default kwargs (modified from [1])
    defaults = (;   formatter=:scientific,
                    legend=:topleft,
                    size=(900, 600),
                    xlabel="Xₚ (km)", 
                    ylabel="Yₚ (km)",
                    title ="Twobody Orbit Positions vs. Time")
    options = merge(defaults, kwargs)

    fig = Plots.plot()

    Plots.plot!(fig, ustrip.(u"km", map(x->perifocal(x)[1], sols)), 
                     ustrip.(u"km", map(x->perifocal(x)[2], sols)), 
                     label="Perifocal Position")
    Plots.plot!(fig; options...)

    return fig

end

function plot3d(sols::Trajectory{<:RestrictedTwoBodyOrbit}; kwargs...)

    # Set default kwargs (modified from [1])
    defaults = (;   formatter=:scientific,
                    legend=:topleft,
                    size=(900, 600),
                    xlabel="X Position (km)", 
                    ylabel="Y Position (km)",
                    zlabel="Z Position (km)",
                    title ="Twobody Orbit Positions vs. Time")
    options = merge(defaults, kwargs)

    fig = Plots.plot()

    Plots.plot!(fig, ustrip.(u"km", map(x->position_vector(x)[1], sols)), 
                     ustrip.(u"km", map(x->position_vector(x)[2], sols)), 
                     ustrip.(u"km", map(x->position_vector(x)[3], sols)), 
                     label="Cartesian Position")
    Plots.plot!(fig; options...)

    return fig

end
