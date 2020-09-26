using Revise, Astrodynamics

r̅₁ = [0.0, 0.0, 0.0]u"km"
v̅₁ = [0.0, 0.0, 0.0]u"km/s"
m₁ = uconvert(u"kg", 1.0u"Mearth")
myEarth = MultibodyState(SVector{3}(r̅₁), SVector{3}(v̅₁), m₁)

r̅₂ = [0.0, 11681.0, 0.0]u"km"
v̅₂ = [5.134, 4.226, 2.787]u"km/s"
m₂ = 150.0u"kg"
mySatellite = MultibodyState(SVector{3}(r̅₂), SVector{3}(v̅₂), m₂)

sys1 = MultibodySystem([myEarth, mySatellite])
sols1 = propagate_multibody(sys1, 10000u"s"; saveat=10.0)

b1 = MultibodyState([ 2.0e6, 0, 0]u"m", [0,  5000.0, 0]u"m/s", 1.0e24u"kg")
b2 = MultibodyState([-2.0e6, 0, 0]u"m", [0, -5000.0, 0]u"m/s", 1.0e24u"kg")
b3 = MultibodyState([ 4.0e6, 0, 0]u"m", [0, -5000.0, 3000.0]u"m/s", 1.0e24u"kg")
b4 = MultibodyState([-4.0e6, 0, 0]u"m", [0,  5000.0, 3000.0]u"m/s", 1.0e24u"kg")

sys2 = MultibodySystem([b1, b2, b3, b4])

sols = propagate_multibody(sys2, 20000.0u"s"; saveat=10.0)

using Plots

plot(map(x-> (x[:,1], x[:,2], x[:,3]), ustrip.(u"km", sols.body[1])))
plot!(map(x-> (x[:,1], x[:,2], x[:,3]), ustrip.(u"km", sols.body[2])))
plot!(title="Multibody Propagation: Position vs. Time",
      label=["Body 1" "Body 2"],
      xlabel="X (km)",
      ylabel="Y (km)",
      zlabel="Z (km)")