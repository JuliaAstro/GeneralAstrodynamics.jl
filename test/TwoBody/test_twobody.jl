module TwoBodyUnitTests

using Test
using UnitfulAstrodynamics

@testset "Transforms" begin
    
    rᵢ = [0.0, 11681.0, 0.0] * u"km"
    vᵢ = [5.134, 4.226, 2.787] * u"km/s"
    orbit = Orbit(rᵢ, vᵢ, Earth) |> KeplerianState

    @test orbit.a == 24509.265399338536 * u"km"
    @test orbit.e == 0.723452708202361
    @test orbit.i == 151.50460766373865 * u"°"
    @test orbit.ν == 89.99652573907436  * u"°"

    @test orbit ≈ TwoBodyState(orbit)

    e      =  0.3
    a      =  15000.   * u"km" + 1.0u"Rearth"
    i      =  10.      * u"°"
    Ω      =  0.       * u"°"
    ω      =  10.      * u"°"
    ν      =  0.       * u"°"
    orbit  =  KeplerianState(e, a, i, Ω, ω, ν, Earth)

    @test isapprox(orbit, TwoBodyState(orbit), atol=1e-6)

end

@testset "Kepler" begin
    
    rᵢ = [0.0, 11681.0, 0.0]u"km"
    vᵢ = [5.134, 4.226, 2.787]u"km/s"
    orbit = Orbit(rᵢ, vᵢ, Earth)

    @test kepler(orbit, period(orbit)) ≈ orbit

end

@testset "Lambert" begin
    
    rᵢ = [0.0, 11681.0, 0.0]u"km"
    vᵢ = [5.134, 4.226, 2.787]u"km/s"
    initial = Orbit(rᵢ, vᵢ, Earth)

    Δt = 1000u"s"
    final = kepler(initial, Δt; tol=1e-12)

    v₁, v₂ = lambert(radius_vector(initial), radius_vector(final), Earth.μ, Δt, :short; tol=1e-12, max_iter=1000)

    @test isapprox(ustrip.(u"km/s", v₁), ustrip.(u"km/s", velocity_vector(initial)); atol=1e-6)
    @test isapprox(ustrip.(u"km/s", v₂), ustrip.(u"km/s", velocity_vector(final)); atol=1e-6)

end

end
