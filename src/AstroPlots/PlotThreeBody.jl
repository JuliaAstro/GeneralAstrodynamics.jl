#
# CR3BP Plots
#

function lagrangeplot(μ, L=1:5; kwargs...)

    defaults = (; title="Nondimensional Lagrange Points", 
                  xlabel="X (DU)", ylabel="Y (DU)", 
                  labels=["Body 1" "Body 2" [string("L",i) for i ∈ L]...])
    options  = merge(defaults, kwargs)

    lagrange_points = lagrange(μ)

    fig = scatter([-μ], [0]; markersize=10, markercolor=:lightblue)
    scatter!(fig, [1-μ], [0]; markersize=6, markercolor=:gray)

    colors = (:red, :orange, :tan, :cyan, :indigo)
    for point ∈ zip(lagrange_points, 1:length(lagrange_points))
        p, i = point
        scatter!(fig, [p[1]], [p[2]]; 
                markershape=:x, markercolor=colors[i], label=string("L", i))
    end

    scatter!(; options...)

    fig

end