using Test
include("TwoBody.jl")

# Test Cartesian -> Canonical transformation
@testset "Transformations" begin
    
    # Set up test initial conditions
    r̅ = [0, 11681, 0] * u"km"
    v̅ = [5.134, 4.226, 2.787] * u"km/s"
    cart = CartesianOrbit(r̅, v̅, earth)

    # Convert Cartesian form to Canonical form
    canon = CanonicalOrbit(cart)

    @test canon.a ≈ 24509.272364065997 * u"km"
    @test canon.e ≈ 0.7234527725236475
    @test canon.i ≈ 2.6442542356744734 * u"rad"
    @test canon.ν ≈ 1.5707355666179315 * u"rad"

    calculated_cart = CartesianOrbit(canon)
    @test calculated_cart ≈ cart

end
