#
# Plot orbits
#

function plot_state(sols::NBodyPropagationResult; kwargs...)

    for i = 1:length(sols.body)

        plot!(ustrip.(u"s", sols.t), ustrip.(u"km", ([sols.body[i].r̅[:,x] for x in 1:size(sols.body[i].r̅,1)]))...)

    end

end