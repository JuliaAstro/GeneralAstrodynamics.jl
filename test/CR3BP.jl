"""
Circular Restricted Three-body Model tests.
"""
module CR3BPTests 

using LinearAlgebra, Unitful, UnitfulAstro, GeneralAstrodynamics, Test

@testset verbose=false "CR3BP Determination" begin
    
    @testset "Unitful" begin
       
        r  = [1.007988, 0.0, 0.001864]u"AU"
        v  = [0, 0, 0]u"AU/s"

        orbit = CR3BPOrbit(r, v, SunEarth)
        @test orbit isa CR3BPOrbit

        orbit = normalize(orbit)
        @test all(position_vector(orbit) .≈ [1.007988, 0.0, 0.001864])
        @test all(velocity_vector(orbit) .≈ zeros(3))

    end

    @testset "No Units" begin
        
        r = [1.2, 0, 0]
        v = [0, -1.049357509830343, 0]
        μ = 0.012150585609624

        orbit = CR3BPOrbit(r, v, μ)

        @test orbit isa NormalizedSynodicCR3BPOrbit

    end

end

@testset verbose=false "CR3BP Calculations" begin
    
    r  = [1.007988, 0.0, 0.001864]u"AU"
    v  = [0, 0, 0]u"AU/s"

    orbit = CR3BPOrbit(r, v, SunEarth) |> normalize

    @test jacobi_constant(orbit) ≈ 3.000907212196274

    @test_broken synodic(inertial(orbit)) ≈ orbit

end

end # module