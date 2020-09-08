using Test
include("TwoBody.jl")

# Test Cartesian -> Canonical transformation
@testset "CartesianToCanonical" begin
    
    # Set up test initial conditions
    r̲ = SVector(0, 11681, 0) * u"km"
    v̲ = SVector(5.134, 4.226, 2.787) * u"km/s"
    cart = CartesianOrbit(r̲, v̲, earth)

    # Convert Cartesian form to Canonical form
    canon = CartesianToCanonical(cart)

    @test canon.a ≈ 24509.272364065997 * u"km"
    @test canon.e ≈ 0.7234527725236475
    @test canon.i ≈ 2.6442542356744734 * u"rad"
    @test canon.ν ≈ 1.5707355666179315 * u"rad"

end
