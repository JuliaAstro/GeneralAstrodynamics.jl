"""
Restricted Two-body Model tests.
"""
module R2BPTests 

using LinearAlgebra, Unitful, GeneralAstrodynamics, Test

@testset verbose=false "R2BP Determination" begin
    
    @testset "Unitful" begin
        rᵢ = [0.0, 11681.0, 0.0] * u"km"
        vᵢ = [5.134, 4.226, 2.787] * u"km/s"
        corbit = CartesianOrbit(rᵢ, vᵢ, Earth, 0.0u"s")

        @test all(
            keplerian(corbit) .≈ (
                0.723452708202361, 
                24509.265399338536u"km", 
                151.50460766373865u"°", 
                90.0u"°", 
                270.0034742609256u"°", 
                89.99652573907436u"°"
            )
        )
        
        @test CartesianOrbit(KeplerianOrbit(corbit)) ≈ corbit

        korbit = KeplerianOrbit(
            0.723452708202361, 
            24509.265399338536u"km", 
            151.50460766373865u"°", 
            90.0u"°", 
            270.0034742609256u"°", 
            89.99652573907436u"°",
            Earth, 
            0.0u"s"
        )

        @test KeplerianOrbit(CartesianOrbit(korbit)) ≈ KeplerianOrbit(corbit)

    end

end

@testset verbose=false "Kepler's Algorithm" begin

    @testset "Unitful" begin
        
        rᵢ = [0.0, 11681.0, 0.0] * u"km"
        vᵢ = [5.134, 4.226, 2.787] * u"km/s"
        corbit = CartesianOrbit(rᵢ, vᵢ, Earth, 0.0u"s")

        @test kepler(corbit) ≈ corbit

    end

end

@testset verbose=true "Lambert Solvers" begin
    
    @testset "Universal" begin
        
        rᵢ = [0.0, 11681.0, 0.0]u"km"
        vᵢ = [5.134, 4.226, 2.787]u"km/s"
        initial = CartesianOrbit(rᵢ, vᵢ, Earth, 0.0u"s")
    
        Δt = 1000u"s"
        final = kepler(initial, Δt; tol=1e-12)
    
        v₁, v₂ = lambert_universal(position_vector(initial), position_vector(final), mass_parameter(Earth), Δt; trajectory=:short, tolerance=1e-6, max_iter=1000)

        @test all(v₁ .≈ velocity_vector(initial))
        @test all(v₂ .≈ velocity_vector(final))

    end

    @testset "Oldenhuis" begin
        
        rᵢ = [0.0, 11681.0, 0.0]u"km"
        vᵢ = [5.134, 4.226, 2.787]u"km/s"
        initial = CartesianOrbit(rᵢ, vᵢ, Earth, 0.0u"s")
    
        Δt = 1000u"s"
        final = kepler(initial, Δt; tol=1e-12)
    
        m = 0
        v₁, v₂ = lambert(position_vector(initial), position_vector(final), Δt, m, mass_parameter(Earth))

        @test_skip all(v₁ .≈ velocity_vector(initial))
        @test_skip all(v₂ .≈ velocity_vector(final))

    end

end


end # module