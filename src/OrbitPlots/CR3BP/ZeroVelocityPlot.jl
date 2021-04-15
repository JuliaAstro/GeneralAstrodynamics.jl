#
# Plot zero velocity curves in the Synodic reference frame with 
# normalized units.
#

"""
Plot the zero velocity curves for the Synodic, normalized CR3BP system.
"""
function zerovelocityplot(orbit::CircularRestrictedThreeBodyOrbit;
                            nondimensional_range = range(-2; stop=2, length=1000),
                            kwargs...)
                                                    
    curves = zerovelocity_curves(orbit; nondimensional_range = nondimensional_range)

    default = (; title  = "Zero Velocity Curves" * (orbit.system.name == "" ? "" : " for the $(orbit.system.name) System"),
                 xlabel = "X ($LU)",
                 ylabel = "Y ($LU)",
                 labels = :none,
                 formatter = :plain,
                 grid      = :on,
                 linewidth = 2)
    options = merge(default, kwargs)

    fig = plot(; options...)
    for curve ∈ curves
        plot!(fig, curve[:,1], curve[:,2])
    end

    return fig
end

"""
Plot the zero velocity curves for the Synodic, normalized CR3BP
system to the last figure.
"""
function zerovelocityplot!(fig, orbit::CircularRestrictedThreeBodyOrbit;
                          nondimensional_range = range(-2; stop=2, length=1000),
                          kwargs...)

    curves = zerovelocity_curves(orbit; nondimensional_range = nondimensional_range)

    default = (; )
    options = merge(default, kwargs)

    fig = plot!(; options...)
    for (i,curve) ∈ zip(1:length(curves), curves)
        plot!(fig, curve[:,1], curve[:,2]; label = "Zero Velocity Curve #$i")
    end

    return fig
end