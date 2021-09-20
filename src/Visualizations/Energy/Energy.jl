#
# Visualizations related to orbital energy
#

# These are not fully implemented!

"""
Plot the zero velocity curves for the Synodic, normalized CR3BP system.
"""
function zerovelocityplot(orbit::CR3BPOrbit;
                          nondimensional_range = range(-2; stop=2, length=1000),
                          kwargs...)
                                                
    curves = zerovelocity_curves(orbit; nondimensional_range = nondimensional_range)

    default = (; title  = "Zero Velocity Curves" * (name(system(orbit)) == (:Primary, :Secondary) ? "" : " for the $(name(system(orbit))[1])-$(name(system(orbit))[2]) System"),
                 xlabel = "X ($(lengthunit(orbit)))",
                 ylabel = "Y ($(lengthunit(orbit)))",
                 labels = :none,
                 formatter = :plain,
                 grid      = :on,
                 linewidth = 2,
                 dpi       = 150)
                                  
    options = merge(default, kwargs)

    fig = plot(; options...)
    for curve ∈ curves
        plot!(fig, curve[:,1], curve[:,2]; kwargs...)
    end

    return fig
end

"""
Plot the zero velocity curves for the Synodic, normalized CR3BP
system to the last figure.
"""
function zerovelocityplot!(fig, orbit::CR3BPOrbit;
                          nondimensional_range = range(-2; stop=2, length=1000),
                          kwargs...)

    curves = zerovelocity_curves(orbit; nondimensional_range = nondimensional_range)

    default = (; )
    options = merge(default, kwargs)

    fig = plot!(; options...)
    for (i,curve) ∈ zip(1:length(curves), curves)
        plot!(fig, curve[:,1], curve[:,2]; kwargs...)
    end

    return fig
end