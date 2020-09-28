#
# Multibody calculations
#

"""
    system_energy(sys::MultibodySystem)

Returns total energy for `NBody` system.
"""
function system_energy(sys::MultibodySystem)

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
    system_angular_momentum(sys::MultibodySystem)

Returns total angular momentum for `NBody` system.
"""
function system_angular_momentum(sys::MultibodySystem)

    H = 0.0u"kg*m^2/s"
    for i = 1:length(sys.body)

        H += sys.body[i].m * cross(sys.body[i].r̅, sys.body[i].v̅)

    end
    return H

end