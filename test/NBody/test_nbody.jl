module NBodyUnitTests

using Test
using UnitfulAstrodynamics

r̅₁ = [0.0, 0.0, 0.0]u"km"
v̅₁ = [0.0, 0.0, 0.0]u"km/s"
m₁ = uconvert(u"kg", 1.0u"Mearth")
myEarth = Body(r̅₁, v̅₁, m₁)

r̅₂ = [0.0, 11681.0, 0.0]u"km"
v̅₂ = [5.134, 4.226, 2.787]u"km/s"
m₂ = 150.0u"kg"
mySatellite = Body(r̅₂, v̅₂, m₂)

sys1 = MultibodySystem([myEarth, mySatellite])

# No errors!
@test true

end
