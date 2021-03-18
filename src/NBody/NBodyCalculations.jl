#
# NBody calculations
#

"""
Returns total energy for `NBody` system.
"""
function system_energy(sys::NBodySystem)

    E = 0.0u"J"
    for i = 1:length(sys.body)

        E += sys.body[i].m * dot(sys.body[i].v̅, sys.body[i].v̅)
        for j = 1:length(sys.body)

            if i ≠ j
                E -= 6.6743e-11u"m^3/(kg*s^2)" * (sys.body[i].m * sys.body[j].m) / 
                                  norm(sys.body[j].r̅ .- sys.body[i].r̅)

            end
        end
    end

    return E

end

"""
Returns total angular momentum for `NBody` system.
"""
function system_angular_momentum(sys::NBodySystem)

    H = reduce(+, [sys.body[i].m * cross(sys.body[i].r̅, sys.body[i].v̅) for i ∈ 1:length(sys.body)])

end