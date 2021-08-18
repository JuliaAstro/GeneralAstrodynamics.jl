#
# Plots for CR3BP manifolds
#

function Plots.plot(manifold::Manifold; vars=:XYZ, kwargs...)

    indices = vars isa Tuple ? vars : process_vars(vars)

    labels = Dict(
        0 => "Time" * " ($(timeunit(manifold)))",
        1 => "X"    * " ($(lengthunit(manifold)))",
        2 => "Y"    * " ($(lengthunit(manifold)))",
        3 => "Z"    * " ($(lengthunit(manifold)))",
        4 => "Ẋ"    * " ($(velocityunit(manifold)))",
        5 => "Ẏ"    * " ($(velocityunit(manifold)))",
        6 => "Ż"    * " ($(velocityunit(manifold)))"
    )

    defaults = (;
        linewidth = 1,
        linestyle = :dot,
        title     = "Manifold",
        palette   = :blues,
        dpi       = 150,
        label     = :none,
        xguide    = labels[indices[1]],
        yguide    = labels[indices[2]],
        zguide    = length(indices) == 3 ? labels[indices[3]] : :none,
    )
    options = merge(defaults, kwargs)
    fig = plot(; options...)

    for trajectory ∈ manifold.solution.u
        plot!(trajectory; vars=vars, options...)
    end

    return fig
end

function Plots.plot!(manifold::Manifold; vars=:XYZ, kwargs...)

    defaults = (;
        linewidth = 1,
        linestyle = :dot,
        palette   = :blues,
        dpi       = 150,
        label     = :none,
    )
    options = merge(defaults, kwargs)

    for trajectory ∈ manifold.solution.u[1:end-1]
        plot!(trajectory; vars=vars, options...)
    end

    plot!(manifold.solution.u[end]; vars=vars, options...)

end

function Plots.plot!(fig, manifold::Manifold; vars=:XYZ, kwargs...)

    defaults = (;
        linewidth = 1,
        linestyle = :dot,
        palette   = :blues,
        dpi       = 150,
        label     = :none,
    )
    options = merge(defaults, kwargs)

    for trajectory ∈ manifold.solution.u
        plot!(fig, trajectory; vars=vars, options...)
    end

    return fig
end