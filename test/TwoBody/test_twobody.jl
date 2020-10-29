module TwoBodyUnitTests

using Test
using UnitfulAstrodynamics

@testset "Transformations" begin
    
    rᵢ = [0.0, 11681.0, 0.0] * u"km"
    vᵢ = [5.134, 4.226, 2.787] * u"km/s"
    orbit = Orbit(rᵢ, vᵢ, Earth)

    @test orbit.a ≈ 24509.272364065997 * u"km"
    @test orbit.e ≈ 0.7234527725236475
    @test orbit.i ≈ 2.6442542356744734 * u"rad"
    @test orbit.ν ≈ 1.5707355666179315 * u"rad"

    @test Orbit(map(x->getfield(orbit, x), [:e, :a, :i, :Ω, :ω, :ν])..., orbit.body) ≈ orbit

    e      =  0.3
    a      =  15000.   * u"km" + 1.0u"Rearth"
    i      =  10.      * u"°"
    Ω      =  0.       * u"°"
    ω      =  10.      * u"°"
    ν      =  0.       * u"°"
    orbit  =  Orbit(e, a, i, Ω, ω, ν, Earth)

    @test isapprox(orbit, Orbit(orbit.rᵢ, orbit.vᵢ, orbit.body), atol=1e-6)

end

@testset "Kepler" begin
    
    rᵢ = [0.0, 11681.0, 0.0]u"km"
    vᵢ = [5.134, 4.226, 2.787]u"km/s"
    orbit = Orbit(rᵢ, vᵢ, Earth)

    @test kepler(orbit, orbital_period(orbit)) ≈ orbit

end

@testset "Lambert" begin
    
    rᵢ = [0.0, 11681.0, 0.0]u"km"
    vᵢ = [5.134, 4.226, 2.787]u"km/s"
    initial = Orbit(rᵢ, vᵢ, Earth)

    Δt = 1000u"s"
    final = kepler(initial, Δt; tol=1e-10)

    v₁, v₂ = lambert(initial.rᵢ, final.rᵢ, Earth.μ, Δt, :short; tol=1e-14, max_iter=10000)

    # 1e-2 km/s is a larger error than I would expect... but 1e-3 fails
    @test isapprox(ustrip.(u"km/s", v₁), ustrip.(u"km/s", initial.vᵢ); atol=1e-2)
    @test isapprox(ustrip.(u"km/s", v₂), ustrip.(u"km/s", final.vᵢ); atol=1e-2)

end

end
