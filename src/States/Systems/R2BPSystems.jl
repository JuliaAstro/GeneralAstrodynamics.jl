# 
# R2BP systems in our solar system
#

"""
Constant `R2BPParameters` for our sun!
"""
const Sun = R2BPParameters(1.327124400419393e11u"km^3/s^2"; name=:Sun)

"""
Constant `R2BPParameters` for Mercury.
"""
const Mercury = R2BPParameters(22031.78000000002u"km^3/s^2"; name=:Mercury)

"""
Constant `R2BPParameters` for Venus.
"""
const Venus = R2BPParameters(324858.592u"km^3/s^2"; name=:Venus)

"""
Constant `R2BPParameters` for your home planet!
"""
const Earth = R2BPParameters(398600.4354360959u"km^3/s^2"; name=:Earth)

"""
Constant `R2BPParameters` for our moon.
"""
const Moon = R2BPParameters(4902.800066163796u"km^3/s^2"; name=:Moon)

"""
Constant `R2BPParameters` (alias for our mooon).
"""
const Luna = Moon

"""
Constant `R2BPParameters` for Mars.
"""
const Mars = R2BPParameters(42828.37362069909u"km^3/s^2"; name=:Mars)

"""
Constant `R2BPParameters` for Jupiter.
"""
const Jupiter = R2BPParameters(1.2668653492180079e8u"km^3/s^2"; name=:Jupiter)

"""
Constant `R2BPParameters` for Saturn.
"""
const Saturn = R2BPParameters(3.793120749865224e7u"km^3/s^2", name=:Saturn)

"""
Constant `R2BPParameters` for Uranus.
"""
const Uranus = R2BPParameters(5.793951322279009e6u"km^3/s^2"; name=:Uranus)

"""
Constant `R2BPParameters` for Neptune.
"""
const Neptune = R2BPParameters(6.835099502439672e6u"km^3/s^2"; name=:Neptune)

"""
Constant `R2BPParameters` for Pluto. We couldn't leave you out again!
"""
const Pluto = R2BPParameters(869.6138177608749u"km^3/s^2"; name=:Pluto)
    