module TwoBodyUnitTests

using Test

using Astrodynamics

@testset "Transformations" begin
    
    r̅ = [0.0, 11681.0, 0.0] * u"km"
    v̅ = [5.134, 4.226, 2.787] * u"km/s"
    cart = CartesianState(r̅, v̅, earth)

    canon = KeplerianState(cart)

    @test canon.a ≈ 24509.272364065997 * u"km"
    @test canon.e ≈ 0.7234527725236475
    @test canon.i ≈ 2.6442542356744734 * u"rad"
    @test canon.ν ≈ 1.5707355666179315 * u"rad"

    @test CartesianState(canon) ≈ cart

    e      =  0.3      * u"rad"
    a      =  15000.   * u"km" + 1.0u"Rearth"
    i      =  10.      * u"°"
    Ω      =  0.       * u"°"
    ω      =  10.      * u"°"
    ν      =  0.       * u"°"
    canon  =  KeplerianState(e, a, i, Ω, ω, ν, earth)

    @test isapprox(canon, KeplerianState(CartesianState(canon)), atol=1e-6)

end

@testset "Propagator" begin
    r̅ = [0.0, 11681, 0.0] * u"km"
    v̅ = [5.134, 4.226, 2.787] * u"km/s"
    cart = CartesianState(r̅, v̅, earth)
    sols = propagate(cart)
    @test CartesianState(sols.r̅[end,:], sols.v̅[end,:], earth) ≈ cart 
end

end