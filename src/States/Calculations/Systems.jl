#
# Methods for orbital systems.
#

"""
Returns the mass parameter of the R2BP system.
"""
massparameter(system::R2BPParameters) = States.get_μ(system) * massparamunit(system)
massparameter(orbit::R2BPOrbit) = massparameter(system(orbit))

"""
Returns the mass parameter of the CR3BP system.
"""
massparameter(system::CR3BPParameters) = States.get_μ(system)
massparameter(orbit::CR3BPOrbit) = massparameter(system(orbit))

"""
Returns the mass parameters of the CR3BP system.
"""
massparameters(system::CR3BPParameters) = (primary_massparameter(system), secondary_massparameter(system))
massparameters(orbit::CR3BPOrbit) = massparameters(system(orbit))

"""
Returns the primary mass parameter of the CR3BP system.
"""
function primary_massparameter(system::CR3BPParameters)
    DU  = lengthunit(system)
    TU  = timeunit(system)
    μ   = massparameter(system)
    ∑μᵢ = DU^3 / ((TU / 2π)^2)

    return ∑μᵢ - μ * ∑μᵢ
end
primary_massparameter(orbit::CR3BPOrbit) = primary_massparameter(system(orbit))

"""
Returns the secondary mass parameter of the CR3BP system.
"""
function secondary_massparameter(system::CR3BPParameters) 
    DU  = lengthunit(system)
    TU  = timeunit(system)
    μ   = massparameter(system)
    ∑μᵢ = DU^3 / ((TU / 2π)^2)

    return μ * ∑μᵢ
end
secondary_massparameter(orbit::CR3BPOrbit) = secondary_massparameter(system(orbit))