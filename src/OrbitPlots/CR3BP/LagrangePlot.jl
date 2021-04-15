#
# Plot the equilibrium solutions for normalized CR3BP dynamics.
#

"""
Plot specified lagrange points in the rotating
reference frame of CR3BP system μ.
"""
function lagrangeplot(μ::Real, L=1:5; kwargs...)

    defaults = (; title  = "Nondimensional Lagrange Points", 
                  xlabel = "X (DU)", ylabel="Y (DU)", 
                  labels = ["Body 1" "Body 2" [string("L",i) for i ∈ L]...],
                  legend = :topleft)
    options  = merge(defaults, kwargs)

    lagrange_points = lagrange(μ)

    fig = scatter([-μ], [0]; markersize=10, markercolor=:lightblue, label="Body 1")
    scatter!(fig, [1-μ], [0]; markersize=6, markercolor=:gray, label="Body 2")

    colors = (:red, :orange, :tan, :cyan, :indigo)
    for point ∈ zip(lagrange_points, 1:length(lagrange_points))
        p, i = point
        scatter!(fig, [p[1]], [p[2]]; 
                markershape=:x, markercolor=colors[i], label=string("L", i))
    end

    scatter!(fig; options...)
    for i ∈ 1:min(length(fig.series_list), length(options.labels))
        fig.series_list[i].plotattributes[:label] = options.labels[i]
    end

    fig

end

"""
Plot specified lagrange points in the rotating
reference frame of CR3BP system μ.
"""
function lagrangeplot(sys::CircularRestrictedThreeBodySystem, L=1:5; normalize = true, exclude_z = true, kwargs...)

    defaults = (; title  = sys.name != "" ? "$(sys.name) Lagrange Points" : "Lagrange Plots", 
                  xlabel = "X ($(string(lengthunit(sys))))", ylabel="Y ($(string(lengthunit(sys))))", 
                  labels = [string("L",i) for i ∈ L], 
                  legend = :topleft)
    options  = merge(defaults, kwargs)

    μ = normalized_mass_parameter(sys)
    lagrange_points = lagrange(μ)

    colors = (:red, :orange, :tan, :cyan, :indigo)
    fig = scatter(; options...)
    for point ∈ zip(lagrange_points, 1:length(lagrange_points))
        p, i = point
        if !normalize 
            p .= ustrip.(lengthunit(sys), p .* normalized_length_unit(sys))
        end
        if exclude_z
            scatter!(fig, [p[1]], [p[2]]; markershape=:x, markercolor=colors[i], markerstrokewidth=2, label=string("L", i))
        else
            scatter!(fig, [p[1]], [p[2], [p[3]]]; markershape=:x, markercolor=colors[i], markerstrokewidth=2, label=string("L", i))
        end
    end

    for i ∈ 1:min(length(fig.series_list), length(options.labels))
        fig.series_list[i].plotattributes[:label] = options.labels[i]
    end

    return fig

end

"""
Plot specified lagrange points in the rotating
reference frame of CR3BP system μ to the last figure.
"""
function lagrangeplot!(fig, sys::CircularRestrictedThreeBodySystem, L=1:5; normalize = true, exclude_z = true, kwargs...)

    defaults = (; )
    options  = merge(defaults, kwargs)

    μ = normalized_mass_parameter(sys)
    lagrange_points = lagrange(μ)

    colors = (:red, :orange, :tan, :cyan, :indigo)
    scatter!(fig; options...)
    for point ∈ zip(lagrange_points, 1:length(lagrange_points))
        p, i = point
        if !normalize 
            p .= ustrip.(lengthunit(sys), p .* normalized_length_unit(sys))
        end
        if exclude_z
            scatter!(fig, [p[1]], [p[2]]; markershape=:x, markercolor=colors[i], markerstrokewidth=2, label=string("L", i))
        else
            scatter!(fig, [p[1]], [p[2], [p[3]]]; markershape=:x, markercolor=colors[i], markerstrokewidth=2, label=string("L", i))
        end
    end

    return fig

end