#
# CR3BP systems in our solar system
#

"""
The Sun-Earth CR3BP system.
"""
const SunVenus = CR3BPParameters(get_μ(Sun) * massparamunit(Sun), get_μ(Venus) * massparamunit(Venus), 108.2e6u"km"; primary=:Sun, secondary=:Venus)

"""
The Sun-Earth CR3BP system.
"""
const SunEarth = CR3BPParameters(get_μ(Sun) * massparamunit(Sun), get_μ(Earth) * massparamunit(Earth), 1.0u"AU"; primary=:Sun, secondary=:Earth)

"""
The Earth-Moon CR3BP system.
"""
const EarthMoon = CR3BPParameters(get_μ(Earth) * massparamunit(Earth), get_μ(Moon) * massparamunit(Moon), 384400u"km"; primary=:Earth, secondary=:Moon)

"""
The Sun-Mars CR3BP system.
"""
const SunMars = CR3BPParameters(get_μ(Sun) * massparamunit(Sun), get_μ(Mars) * massparamunit(Mars), 227.9e6u"km"; primary=:Sun, secondary=:Mars)

"""
The Sun-Jupiter CR3BP system.
"""
const SunJupiter = CR3BPParameters(get_μ(Sun) * massparamunit(Sun), get_μ(Jupiter) * massparamunit(Jupiter), 778.6e6u"km"; primary=:Sun, secondary=:Jupiter)

"""
The Sun-Saturn CR3BP system.
"""
const SunSaturn = CR3BPParameters(get_μ(Sun) * massparamunit(Sun), get_μ(Saturn) * massparamunit(Saturn), 1433.5e6u"km"; primary=:Sun, secondary=:Saturn)
