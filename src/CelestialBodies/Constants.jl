#
# Solar System Constants
# 



# Constants

# All data pulled from the following references:
# [1] https://en.wikipedia.org/wiki/List_of_Solar_System_objects_by_size
# [2] https://docs.astropy.org/en/stable/constants/#module-astropy.constants

"""
Constant `CelestialBody` for our sun!
"""
const Sun = CelestialBody(1.327124400419393e11u"km^3/s^2", 696000.0u"km", "Sun")

"""
Constant `CelestialBody` for Mercury.
"""
const Mercury = CelestialBody(22031.78000000002u"km^3/s^2", 2439.7u"km", "Mercury")

"""
Constant `CelestialBody` for Venus.
"""
const Venus = CelestialBody(324858.592u"km^3/s^2", 6051.8u"km", "Venus")

"""
Constant `CelestialBody` for your home planet!
"""
const Earth = CelestialBody(398600.4354360959u"km^3/s^2", 6371.008366666666u"km", "Earth")

"""
Constant `CelestialBody` for our moon.
"""
const Moon = CelestialBody(4902.800066163796u"km^3/s^2", 1737.4000000000003u"km", "Moon")

"""
Constant `CelestialBody` (alias for our mooon).
"""
const Luna = Moon

"""
Constant `CelestialBody` for Mars.
"""
const Mars = CelestialBody(42828.37362069909u"km^3/s^2", 3389.5266666666666u"km", "Mars")

"""
Constant `CelestialBody` for Jupiter.
"""
const Jupiter = CelestialBody(1.2668653492180079e8u"km^3/s^2", 69946.0u"km", "Jupiter")

"""
Constant `CelestialBody` for Saturn.
"""
const Saturn = CelestialBody(3.793120749865224e7u"km^3/s^2", 58300.0u"km", "Saturn")

"""
Constant `CelestialBody` for Uranus.
"""
const Uranus = CelestialBody(5.793951322279009e6u"km^3/s^2", 25363.666666666668u"km", "Uranus")

"""
Constant `CelestialBody` for Neptune.
"""
const Neptune = CelestialBody(6.835099502439672e6u"km^3/s^2", 24623.0u"km", "Neptune")

"""
Constant `CelestialBody` for Pluto. We couldn't leave you out again!
"""
const Pluto = CelestialBody(869.6138177608749u"km^3/s^2", 1195.0u"km", "Pluto")
    
