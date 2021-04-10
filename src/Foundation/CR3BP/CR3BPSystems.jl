#
# Circular Restricted Three-body Problem Systems
# 

"""
The Sun-Earth CR3BP system.
"""
const SunEarth = CircularRestrictedThreeBodySystem(mass_parameter.((Sun, Earth)), 1.0u"AU", time_scale_factor(1.0u"AU", mass_parameter.((Sun, Earth))...), "Sun-Earth")
