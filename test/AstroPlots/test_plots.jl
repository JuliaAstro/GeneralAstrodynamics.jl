module PlotsUnitTests

using Test
using UnitfulAstrodynamics

@testset "AstroPlots" begin
    
    # Twobody Orbit
    r̅ = [0.0, 11681, 0.0] * u"km"
    v̅ = [5.134, 4.226, 2.787] * u"km/s"
    orbit = Orbit(r̅, v̅, Earth)

    # Does plot run?
    fig = plot(propagate(orbit, 1u"s"))

    @test true

end

end
