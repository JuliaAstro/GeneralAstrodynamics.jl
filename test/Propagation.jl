"""
Circular Restricted Three-body Model tests.
"""
module PropagationTests 

using LinearAlgebra, Unitful, UnitfulAstro, GeneralAstrodynamics, Test

@testset verbose=false "R2BP Propagation" begin
    
    orbit = KeplerianOrbit(0.4, 10_000, 0, 0, 0, 0, Mars, 0.0) |> CartesianOrbit
    @test propagate(orbit) ≈ orbit

end

@testset verbose=false "CR3BP Propagation" begin
    
    r  = [1.007988, 0.0, 0.001864]u"AU"
    v  = [0, 0, 0]u"AU/s"

    orbit = CR3BPOrbit(r, v, SunEarth) |> normalize

    T = 3.0967384408950527
    final = propagate(orbit, T)[end]
    @test_broken final ≈ orbit

    @test all(position_vector(final) .≈ position_vector(orbit))
    @test all(Velocity_vector(final) .≈ velocity_vector(orbit))

end

@testset verbose=false "Halo Solvers" begin
    
    @test isperiodic(halo(SunEarth;   L=2, Az=0.0005)...)
    @test isperiodic(halo(SunMars;    L=1, Az=0.0005)...)
    @test isperiodic(halo(EarthMoon;  L=2, Az=0.0005)...)
    @test isperiodic(halo(SunJupiter; L=1, Az=0.0005)...)
    
    orbit, T = halo(SunJupiter; L=2, Az=0)
    @test monodromy(orbit, T) ≈ [
         1177.62   -43.3137  0.0         247.085   227.411   0.0
        -1086.0     40.8647  0.0        -227.411  -209.976   0.0
            0.0      0.0     0.996642      0.0       0.0    -0.0918265
        3446.16  -126.529   0.0         722.795   666.042   0.0
        -2146.97    79.0345  0.0        -450.857  -413.957   0.0
            0.0      0.0     0.0730094     0.0       0.0     0.996642
    ]
    
end

end # module