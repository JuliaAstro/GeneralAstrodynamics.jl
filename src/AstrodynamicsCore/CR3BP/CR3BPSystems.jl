#
# Circular Restricted Three-body Problem Systems
# 

"""
The Sun-Earth CR3BP system.
"""
const SunVenus = CircularRestrictedThreeBodySystem(
    mass_parameter.((Sun, Venus)), 
    108.2e6u"km", 
    time_scale_factor(108.2e6u"km", mass_parameter.((Sun, Venus))...), 
    "Sun-Venus")

"""
The Sun-Earth CR3BP system.
"""
const SunEarth = CircularRestrictedThreeBodySystem(
    mass_parameter.((Sun, Earth)), 
    1.0u"AU", 
    time_scale_factor(1.0u"AU", mass_parameter.((Sun, Earth))...), 
    "Sun-Earth")

"""
The Earth-Moon CR3BP system.
"""
const EarthMoon = CircularRestrictedThreeBodySystem(
    mass_parameter.((Earth, Moon)), 
    384400u"km", 
    time_scale_factor(384400u"km", mass_parameter.((Earth, Moon))...), 
    "Earth-Moon")

"""
The Sun-Mars CR3BP system.
"""
const SunMars = CircularRestrictedThreeBodySystem(
    mass_parameter.((Sun, Mars)), 
    227.9e6u"km", 
    time_scale_factor(227.9e6u"km", mass_parameter.((Sun, Mars))...), 
    "Sun-Mars")

"""
The Sun-Jupiter CR3BP system.
"""
const SunJupiter = CircularRestrictedThreeBodySystem(
    mass_parameter.((Sun, Jupiter)), 
    778.6e6u"km", 
    time_scale_factor(778.6e6u"km", mass_parameter.((Sun, Jupiter))...), 
    "Sun-Jupiter")

"""
The Sun-Saturn CR3BP system.
"""
const SunSaturn = CircularRestrictedThreeBodySystem(
    mass_parameter.((Sun, Saturn)), 
    1433.5e6u"km", 
    time_scale_factor(1433.5e6u"km", mass_parameter.((Sun, Saturn))...), 
    "Sun-Saturn")

