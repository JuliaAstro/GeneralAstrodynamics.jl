#
# Multibody calculations
#

function system_energy(sys::MultibodySystem)

    E = 0.0u"J"

    for i = 1:length(sys.body)

        E += sys.body[i] * dot(sys.v̅, sys.v̅)
        for j = 1:length(sys.body)

            if i ≠ j
                E -= 6.6743e-11 * (sys.body[i].m * sys.body[j].m) / norm(sys.body[j].r̅ .- sys.body[i].r̅)

            end
        end

    end

    return E
    
end

function system_angular_momentum(sys::MultibodySystem)

    H = 0.0u"kg*m/s"

    for i = 1:length(sys.body)

        H += sys.body[i].m * cross(sys.body[i].r̅, sys.body[i].v̅)

    end

    return H

end