#
# Plot orbits
#
# Referencing:
# [1] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

"""
    plot(sols::TwobodyPropagationResult, frame=:Perifocal; kwargs...)

Plots every timestep in `sols` in `3D` space. All keyward 
arguments are passed directly to `Plots.jl`.
"""
function plot(sols::TwobodyPropagationResult, frame=:Perifocal; kwargs...)
   
    # Provided frame can be :Cartesian, or :Perifocal
    if frame == :Perifocal
        return plot2d(sols; kwargs...)
    elseif frame == :Cartesian
        return plot3d(sols; kwargs...)
    else 
        throw(ArgumentError("`frame` must be set to `:Cartesian` or `:Perifocal`"))
    end

end

function plot2d(sols::TwobodyPropagationResult; kwargs...)

    # Set default kwargs (modified from [1])
    defaults = (;   formatter=:scientific,
                    legend=:topleft,
                    size=(900, 600),
                    xlabel="Xₚ (km)", 
                    ylabel="Yₚ (km)",
                    title ="Twobody Orbit Positions vs. Time")
    options = merge(defaults, kwargs)

    fig = Plots.plot()

    Plots.plot!(fig, ustrip.(u"km", map(x->x.rₚ[1], sols.step)), 
                     ustrip.(u"km", map(x->x.rₚ[2], sols.step)), 
                     label="Perifocal Position")
    Plots.plot!(fig; options...)

    return fig

end

function plot3d(sols::TwobodyPropagationResult; kwargs...)

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

    Plots.plot!(fig, ustrip.(u"km", map(x->x.rᵢ[1], sols.step)), 
                     ustrip.(u"km", map(x->x.rᵢ[2], sols.step)), 
                     ustrip.(u"km", map(x->x.rᵢ[3], sols.step)), 
                     label="Cartesian Position")
    Plots.plot!(fig; options...)

    return fig

end
