#
# Defines AstrodynamicalCalculations methods for 
# OrbitalStates types
#

export massparameter, massparameters, normalized_massparameter
export primary_massparameter, secondary_massparameter
export RAAN, argument_of_periapsis, inclination
export periapsis_velocity, apoapsis_velocity
export mean_motion_vector
export distance_to_primary, distance_to_secondary
export massparameter, massparameters
export primary_massparameter, secondary_massparameter

include("States.jl")
include("Systems.jl")