"""
Circular Restricted Three-body Model tests.
"""
module CR3BPTests 

using LinearAlgebra, Unitful, UnitfulAstro, GeneralAstrodynamics, Test

@testset verbose=false "CR3BP Determination" begin
    
    @testset "Unitful" begin
       
        r  = [1.007988, 0.0, 0.001864]u"AU"
        v  = [0, 0, 0]u"AU/s"

        orbit = Orbit(CartesianState(r, v), SunEarth)

        @test all(orbit.state.r .≈ [1.007988, 0.0, 0.001864])
        @test all(orbit.state.v .≈ zeros(3))

    end

    @testset "No Units" begin
        
        r = [1.2, 0, 0]
        v = [0, -1.049357509830343, 0]
        μ = 0.012150585609624

        units = (; lengthunit=missing, timeunit=missing, angularunit=missing)
        @test_broken Orbit(
            CartesianState(r, v; units...), 
            CR3BPParameters(μ; massunit=missing, units...)
        ).state ≈ vcat(r,v)

    end

end

@testset verbose=false "CR3BP Calculations" begin
    
    r  = [1.007988, 0.0, 0.001864]u"AU"
    v  = [0, 0, 0]u"AU/s"

    orbit = Orbit(CartesianState(r, v), SunEarth) 

    @test jacobi_constant(orbit) ≈ 3.000907212196274

    @test_broken synodic(inertial(orbit)) ≈ orbit

end

end # module