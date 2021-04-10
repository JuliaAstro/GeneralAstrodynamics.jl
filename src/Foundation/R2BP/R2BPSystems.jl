#
# Provides out-of-the-box Two-body systems in our solar system.
#

"""
Constant `RestrictedTwoBodySystem` for our sun!
"""
const Sun = RestrictedTwoBodySystem(1.327124400419393e11u"km^3/s^2", "Sun")

"""
Constant `RestrictedTwoBodySystem` for Mercury.
"""
const Mercury = RestrictedTwoBodySystem(22031.78000000002u"km^3/s^2", "Mercury")

"""
Constant `RestrictedTwoBodySystem` for Venus.
"""
const Venus = RestrictedTwoBodySystem(324858.592u"km^3/s^2", "Venus")

"""
Constant `RestrictedTwoBodySystem` for your home planet!
"""
const Earth = RestrictedTwoBodySystem(398600.4354360959u"km^3/s^2", "Earth")

"""
Constant `RestrictedTwoBodySystem` for our moon.
"""
const Moon = RestrictedTwoBodySystem(4902.800066163796u"km^3/s^2", "Moon")

"""
Constant `RestrictedTwoBodySystem` (alias for our mooon).
"""
const Luna = Moon

"""
Constant `RestrictedTwoBodySystem` for Mars.
"""
const Mars = RestrictedTwoBodySystem(42828.37362069909u"km^3/s^2", "Mars")

"""
Constant `RestrictedTwoBodySystem` for Jupiter.
"""
const Jupiter = RestrictedTwoBodySystem(1.2668653492180079e8u"km^3/s^2", "Jupiter")

"""
Constant `RestrictedTwoBodySystem` for Saturn.
"""
const Saturn = RestrictedTwoBodySystem(3.793120749865224e7u"km^3/s^2", "Saturn")

"""
Constant `RestrictedTwoBodySystem` for Uranus.
"""
const Uranus = RestrictedTwoBodySystem(5.793951322279009e6u"km^3/s^2", "Uranus")

"""
Constant `RestrictedTwoBodySystem` for Neptune.
"""
const Neptune = RestrictedTwoBodySystem(6.835099502439672e6u"km^3/s^2", "Neptune")

"""
Constant `RestrictedTwoBodySystem` for Pluto. We couldn't leave you out again!
"""
const Pluto = RestrictedTwoBodySystem(869.6138177608749u"km^3/s^2", "Pluto")
    