module PropagatorUnitTests

using Test
using UnitfulAstrodynamics

@testset "Propagators" begin

    # Twobody Orbit
    r̅ = [0.0, 11681, 0.0] * u"km"
    v̅ = [5.134, 4.226, 2.787] * u"km/s"
    orbit = Orbit(r̅, v̅, Earth)

    # Propagate twobody
    sols = propagate(orbit, 5u"s"; save_everystep=false)
    @test isapprox(sols.step[end], kepler(orbit, 5u"s"; tol=1e-14), atol=1e-6) 

end

end
