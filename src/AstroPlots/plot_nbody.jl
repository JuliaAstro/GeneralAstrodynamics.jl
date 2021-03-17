#
# Plot orbits
#

"""
Plots every timestep in `sols` in `3D` space. All keyward 
arguments are passed directly to `Plots.jl`.
"""
function orbitplot(sols::MultibodyPropagationResult; bodies=1:length(sols.step[1].body), kwargs...)
   
    # Referencing:
    # [1] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

    # Set default kwargs (modified from [1])
    defaults = (;   formatter=:scientific,
                    legend=:topleft,
                    xlabel="X Position (km)", 
                    ylabel="Y Position (km)",
                    zlabel="Z Position (km)",
                    title ="NBody Positions vs. Time")
    options = merge(defaults, kwargs)

    fig = Plots.plot()
    for i = bodies

        Plots.plot!(fig, ustrip.(u"km", map(x -> x.body[i].r̅[1], sols.step)), 
                         ustrip.(u"km", map(x -> x.body[i].r̅[2], sols.step)), 
                         ustrip.(u"km", map(x -> x.body[i].r̅[3], sols.step)), 
                         label=string("Body ", i))

    end
    Plots.plot!(fig; options...)

    return fig
    
end