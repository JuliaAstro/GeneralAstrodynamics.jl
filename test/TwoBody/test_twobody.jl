module TwoBodyUnitTests

using Test
using UnitfulAstrodynamics
using UnitfulAstrodynamics.TwoBody.Systems: Earth

@testset "Transforms" begin
    
    rᵢ = [0.0, 11681.0, 0.0] * u"km"
    vᵢ = [5.134, 4.226, 2.787] * u"km/s"
    orbit = CartesianOrbit(rᵢ, vᵢ, Earth, 0.0u"s") |> KeplerianOrbit

    @test semimajor_axis(orbit) == 24509.265399338536 * u"km"
    @test eccentricity(orbit) == 0.723452708202361
    @test inclination(orbit) == 151.50460766373865 * u"°"
    @test true_anomoly(orbit) == 89.99652573907436  * u"°"

    @test orbit ≈ CartesianOrbit(orbit)

    e      =  0.3
    a      =  semimajor_axis(orbit)
    i      =  10.      * u"°"
    Ω      =  0.       * u"°"
    ω      =  10.      * u"°"
    ν      =  0.       * u"°"
    orbit  =  KeplerianOrbit(e, a, i, Ω, ω, ν, Earth, 0.0u"s")

    @test isapprox(orbit, CartesianOrbit(orbit); atol=1e-6)

end

@testset "Kepler" begin
    
    rᵢ = [0.0, 11681.0, 0.0]u"km"
    vᵢ = [5.134, 4.226, 2.787]u"km/s"
    orbit = CartesianOrbit(rᵢ, vᵢ, Earth, 0.0u"s")

    @test kepler(orbit) ≈ orbit

end

@testset "Lambert" begin
    
    rᵢ = [0.0, 11681.0, 0.0]u"km"
    vᵢ = [5.134, 4.226, 2.787]u"km/s"
    initial = CartesianOrbit(rᵢ, vᵢ, Earth, 0.0u"s")

    Δt = 1000u"s"
    final = kepler(initial, Δt; tol=1e-12)

    v₁, v₂ = lambert(position_vector(initial), position_vector(final), mass_parameter(Earth), Δt, :short; tol=1e-6, max_iter=1000)

    @test isapprox(ustrip.(u"km/s", v₁), ustrip.(u"km/s", velocity_vector(initial)); atol=1e-6)
    @test isapprox(ustrip.(u"km/s", v₂), ustrip.(u"km/s", velocity_vector(final)); atol=1e-6)

end

end
