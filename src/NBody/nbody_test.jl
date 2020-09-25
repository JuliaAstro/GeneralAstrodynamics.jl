using Revise, Astrodynamics, ComponentArrays

r̅₁ = [0.0, 0.0, 0.0]u"km"
v̅₁ = [0.0, 0.0, 0.0]u"km/s"
m₁ = uconvert(u"kg", 1.0u"Mearth")
myEarth = BodyState(SVector{3}(r̅₁), SVector{3}(v̅₁), m₁)

r̅₂ = [0.0, 11681.0, 0.0]u"km"
v̅₂ = [5.134, 4.226, 2.787]u"km/s"
m₂ = 150.0u"kg"
mySatellite = BodyState(SVector{3}(r̅₂), SVector{3}(v̅₂), m₂)

sys1 = System([myEarth, mySatellite])

B = Array([ComponentArray((r̅=r̅₁, v̅=v̅₁)), ComponentArray((r̅=r̅₂, v̅=v̅₂))])
p = ComponentArray((G=6.6743e-11, m=[m₁, m₂]))

b1 = BodyState([2.0e6, 0, 0]u"m", [0, 5000.0, 0]u"m/s", 1.0e24u"kg")
b2 = BodyState([-2.0e6, 0, 0]u"m", [0, -5000.0, 0]u"m/s", 1.0e24u"kg")
b3 = BodyState([4.0e6, 0, 0]u"m", [0, -5000.0, 3000.0]u"m/s", 1.0e24u"kg")
b4 = BodyState([-4.0e6, 0, 0]u"m", [0, 5000.0, 3000.0]u"m/s", 1.0e24u"kg")

sys2 = System([b1, b2, b3, b4])

sols = npropagate(sys2, 20000.0u"s"; saveat=10.0)