module PlotsUnitTests

using Test
using UnitfulAstrodynamics
using UnitfulAstrodynamics.TwoBody.Systems: Earth
@testset "OrbitPlots" begin
    
    # Twobody Orbit
    r̅ = [0.0, 11681, 0.0] * u"km"
    v̅ = [5.134, 4.226, 2.787] * u"km/s"
    orbit = CartesianOrbit(r̅, v̅, Earth, 0.0u"s")

    # Does plot run?
    fig = orbitplot(propagate(orbit, 1u"s"))

    @test true

end

end
